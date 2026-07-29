#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
})

parse_args <- function(x) {
  out <- list(
    input = NULL,
    gtf = NULL,
    outdir = "gene_variant_plots",
    width = 10,
    height = 6,
    dpi = 300
  )
  i <- 1
  while (i <= length(x)) {
    flag <- x[[i]]
    val <- if (i < length(x)) x[[i + 1]] else NA_character_
    if (flag %in% c("-i", "--input")) {
      out$input <- val
      i <- i + 2
    } else if (flag %in% c("-g", "--gtf")) {
      out$gtf <- val
      i <- i + 2
    } else if (flag %in% c("-o", "--outdir")) {
      out$outdir <- val
      i <- i + 2
    } else if (flag == "--width") {
      out$width <- as.numeric(val)
      i <- i + 2
    } else if (flag == "--height") {
      out$height <- as.numeric(val)
      i <- i + 2
    } else if (flag == "--dpi") {
      out$dpi <- as.numeric(val)
      i <- i + 2
    } else if (flag %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  Rscript 8_plot_gene_variants.R --input summary.tsv --gtf genes.gtf --outdir plots\n\n",
        "Required input columns:\n",
        "  gene_id, gene_name, variant_ids, genotypes, participant_id, case, dataset\n\n",
        "Notes:\n",
        "  - variant_ids and genotypes may contain semicolon-separated values on the same row\n",
        "  - the GTF is used to look up gene start/end from gene features\n",
        sep = ""
      )
      quit(status = 0)
    } else {
      stop(sprintf("Unknown argument: %s", flag))
    }
  }
  if (is.null(out$input)) stop("Missing required --input")
  out
}

split_semicolon <- function(x) {
  parts <- str_trim(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  parts[nzchar(parts)]
}

normalize_gt <- function(gt) {
  str_replace_all(str_trim(as.character(gt)), "\\|", "/")
}

parse_variant_pos <- function(variant_id) {
  parts <- strsplit(as.character(variant_id), ":", fixed = TRUE)[[1]]
  if (length(parts) != 4) stop(sprintf("Variant id must look like chr:pos:ref:alt, got: %s", variant_id))
  as.integer(parts[[2]])
}

safe_name <- function(x) {
  x <- str_replace_all(as.character(x), "[^A-Za-z0-9_.-]+", "_")
  str_replace_all(x, "_+", "_")
}

parse_gtf_attributes <- function(attrs) {
  fields <- unlist(strsplit(attrs, ";", fixed = TRUE))
  fields <- str_trim(fields)
  fields <- fields[nzchar(fields)]
  out <- list()
  for (field in fields) {
    if (!str_detect(field, " ")) next
    key <- str_trim(str_split_fixed(field, " ", 2)[1, 1])
    value <- str_trim(str_split_fixed(field, " ", 2)[1, 2])
    out[[key]] <- str_remove_all(value, '"')
  }
  out
}

read_gtf_gene_coordinates <- function(gtf_path) {
  con <- file(gtf_path, open = "r")
  on.exit(close(con), add = TRUE)
  rows <- list()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (startsWith(line, "#")) next
    parts <- str_split_fixed(line, "\t", 9)
    if (parts[1, 3] != "gene") next
    attrs <- parse_gtf_attributes(parts[1, 9])
    if (is.null(attrs$gene_id)) next
    rows[[length(rows) + 1]] <- tibble(
      gene_id = attrs$gene_id,
      gene_name = if (!is.null(attrs$gene_name)) attrs$gene_name else attrs$gene_id,
      gene_start = as.integer(parts[1, 4]),
      gene_end = as.integer(parts[1, 5]),
      strand = parts[1, 7]
    )
  }
  bind_rows(rows) %>% distinct(gene_id, .keep_all = TRUE)
}

expand_rows <- function(df) {
  map_dfr(seq_len(nrow(df)), function(i) {
    row <- df[i, , drop = FALSE]
    vids <- split_semicolon(row$variant_ids[[1]])
    gts <- split_semicolon(row$genotypes[[1]])
    if (length(vids) != length(gts)) {
      stop(sprintf(
        "Row %d has %d variant_ids but %d genotypes",
        i, length(vids), length(gts)
      ))
    }
    tibble(
      gene_id = row$gene_id[[1]],
      gene_name = row$gene_name[[1]],
      participant_id = row$participant_id[[1]],
      case = row$case[[1]],
      dataset = row$dataset[[1]],
      variant_id = vids,
      genotype = gts
    )
  })
}

assign_gene_coordinates <- function(dat, gtf_path) {
  if (is.null(gtf_path)) {
    stop("Missing required --gtf")
  }
  ann <- read_gtf_gene_coordinates(gtf_path)
  joined <- left_join(dat, ann, by = "gene_id", suffix = c("", "_gtf"))
  if (all(is.na(joined$gene_start))) {
    stop("GTF join did not provide gene_start for any rows")
  }
  joined
}

make_gene_plot <- function(dat, out_path) {
  dat <- dat %>%
    mutate(
      variant_pos = map_int(variant_id, parse_variant_pos),
      genotype_norm = normalize_gt(genotype),
      gt_class = case_when(
        genotype_norm %in% c("0/1", "1/0") ~ "het",
        genotype_norm == "1/1" ~ "hom",
        TRUE ~ "other"
      ),
      x_pos = variant_pos
    )

  participant_levels <- dat %>%
    distinct(participant_id, case) %>%
    mutate(case_rank = if_else(case == "control", 1L, 2L, missing = 3L)) %>%
    arrange(case_rank, participant_id) %>%
    pull(participant_id)

  dat <- dat %>%
    mutate(participant_id = factor(participant_id, levels = participant_levels))

  segments <- dat %>%
    group_by(gene_id, gene_name, dataset, participant_id, case) %>%
    summarise(
      xmin = min(x_pos, na.rm = TRUE),
      xmax = max(x_pos, na.rm = TRUE),
      .groups = "drop"
    )

  segments <- segments %>%
    mutate(participant_id = factor(participant_id, levels = participant_levels))

  plot_title <- sprintf("%s (%s)", first(dat$gene_name), first(dat$gene_id))
  plot_subtitle <- paste0("Dataset: ", paste(unique(dat$dataset), collapse = ", "))

  p <- ggplot() +
    geom_segment(
      data = segments,
      aes(x = xmin, xend = xmax, y = participant_id, yend = participant_id, colour = case),
      linewidth = 1
    ) +
    geom_point(
      data = filter(dat, gt_class == "het"),
      aes(x = x_pos, y = participant_id, colour = case),
      shape = 1,
      size = 2.8,
      stroke = 0.9
    ) +
    geom_point(
      data = filter(dat, gt_class == "hom"),
      aes(x = x_pos, y = participant_id, colour = case),
      shape = 16,
      size = 2.8
    ) +
    geom_point(
      data = filter(dat, gt_class == "other"),
      aes(x = x_pos, y = participant_id, colour = case),
      shape = 4,
      size = 2.8,
      stroke = 0.9
    ) +
    scale_colour_manual(values = c(case = "#D1495B", control = "#2E86AB")) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Genomic position (bp)",
      y = "Participant",
      colour = "Group"
    ) +
    scale_x_continuous(
      limits = c(min(dat$gene_start, na.rm = TRUE), max(dat$gene_end, na.rm = TRUE)),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave(out_path, p, width = 10, height = max(4, 0.25 * length(participant_levels) + 2), dpi = 300)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read_tsv(args$input, show_col_types = FALSE)
required <- c("gene_id", "gene_name", "variant_ids", "genotypes", "participant_id", "case", "dataset")
missing <- setdiff(required, names(dat))
if (length(missing) > 0) {
  stop(sprintf("Input TSV is missing columns: %s", paste(missing, collapse = ", ")))
}

dat <- expand_rows(dat)
dat <- assign_gene_coordinates(dat, args$gtf)

plots <- dat %>%
  group_by(gene_id, gene_name) %>%
  group_split()

walk(plots, function(df) {
  gene_id <- as.character(first(df$gene_id))
  gene_name <- as.character(first(df$gene_name))
  file_base <- safe_name(sprintf("%s__%s", gene_id, gene_name))
  out_path <- file.path(args$outdir, paste0(file_base, ".pdf"))
  message("Writing ", out_path)
  make_gene_plot(df, out_path)
})
