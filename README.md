# snoRNA_biallelic

Minimal workflow for finding individuals with multiple snoRNA variants in AGGV3 VCF shards.

## What this repo does

- Uses a GTF to define snoRNA intervals for coverage scoring.
- Calculates a coverage score from AGGV3 site-QC VCFs before variant querying.
- Uses the coverage summary as the gene input for variant querying, so step 2 does not reopen the GTF.
- Uses `biallelic_shards.bed` to find the relevant shard/subshard VCFs for each snoRNA gene.
- `2_write_gene_variant_tsvs.py` skips genes below the coverage threshold and writes one TSV per gene as soon as that gene finishes.
- Writes rare variant rows only, with participant genotype plus `AF`, `AC`, and `AN` in each per-gene TSV.
- Merges nearby genes into shared fetch windows so each shard is queried fewer times.

## Suggested layout on CloudOS

- `annotations.gtf.gz`
- `vcfs/`
- `biallelic_shards.bed`
- `participants.tsv` optional, if you want downstream participant annotation

## Coverage First

To calculate a coverage score per snoRNA from the AGGV3 site-QC VCFs, run:

```bash
python scripts/1_calculate_snoRNA_coverage_score.py \
  --gtf annotations.gtf.gz \
  --shard-bed biallelic_shards.bed \
  --site-qc-root site_qc_vcfs \
  --vcf-template shard-{shard}/subshard-{subshard}/postproc/site_qc.vcf.gz \
  --out outputs/snorna_biallelic.coverage_score.tsv
```

This uses cohort `MEDIAN_DP` from the site-QC VCF records that overlap each snoRNA interval, so it is a proxy coverage score rather than a true base-wise depth track.

Interpretation:

- `n_site_qc_records` is the number of site-QC records overlapping the snoRNA that had a usable `MEDIAN_DP` value.
- `coverage_score` is the mean `MEDIAN_DP` across those overlapping records.
- If `n_site_qc_records` is `0`, then `coverage_score` will be `0.0`.

## Variant Query

```bash
python scripts/2_write_gene_variant_tsvs.py \
  --coverage-summary outputs/snorna_biallelic.coverage_score.tsv \
  --min-coverage-score 20 \
  --shard-bed biallelic_shards.bed \
  --vcf-root vcfs \
  --cpus 16 \
  --out-prefix outputs/snorna_biallelic
```

This writes:

- `outputs/snorna_biallelic.genes/0001_*.tsv`

To find participants with rare homozygous snoRNA variants, run:

```bash
python scripts/3_find_hom_vars.py \
  --genes-dir outputs/snorna_biallelic.genes \
  --out outputs/snorna_biallelic.homozygous_variants.tsv
```

To find participants with at least two rare heterozygous variants in the same snoRNA, run:

```bash
python scripts/4_find_double_het_vars.py \
  --genes-dir outputs/snorna_biallelic.genes \
  --out outputs/snorna_biallelic.two_rare_same_snoRNA.tsv
```

To liftover GRCh38 variant IDs to GRCh37 while preserving the input TSV and adding a lifted variant column, run:

```bash
python scripts/5_liftover_variant_ids_grch38_to_grch37.py \
  --variants variants.grch38.tsv \
  --column variant_id \
  --chain hg38ToHg19.over.chain.gz \
  --fasta GRCh37.fa \
  --out variants.grch37.tsv \
  --unmapped variants.unmapped.tsv
```

This expects the UCSC `liftOver` executable on `PATH` and a GRCh37 FASTA with a matching `.fai` index.

If you run script 2 and your mounted directory structure differs from the default `shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz` pattern, pass `--vcf-template` with the relative path layout that matches your session.

Script 2 prints progress as it loads genes, queues each shard, and reports each shard as it finishes.
It prefers `cyvcf2` for indexed region fetches when available, then falls back to `tabix -h`, and finally to a plain scan of the VCF file. If you know the VCFs are indexed and `cyvcf2` is installed, `--region-access auto` is the fastest option.

## Participant TSV

The participant table is not needed by scripts 1 to 4, but you can keep one around for later case/control annotation or other downstream filtering.

## Notes

- The script runs without extra dependencies, but `cyvcf2` is used automatically when available.
- It treats a participant as a carrier for a site if any genotype allele is non-reference.
- Multi-allelic VCF rows are counted once per gene and participant.
- Final participant and gene summaries include only variants with `AF < 0.005`.
