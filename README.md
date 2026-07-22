# snoRNA_biallelic

Minimal workflow for finding individuals with multiple snoRNA variants in AGGV3 VCF shards.

## What this repo does

- Reads a GTF and builds snoRNA intervals from `gene_type` or `gene_biotype`.
- Uses `biallelic_shards.bed` to find the relevant shard/subshard VCFs for each snoRNA gene.
- `1_count_snoRNA_multi_variant_carriers.py` writes one TSV per gene as soon as that gene finishes.
- Produces end-of-run gene and participant summaries for rare variants only.
- Records the participant genotype plus `AF`, `AC`, and `AN` in each per-gene TSV.
- Merges nearby genes into shared fetch windows so each shard is queried fewer times.

## Suggested layout on CloudOS

- `annotations.gtf.gz`
- `vcfs/`
- `biallelic_shards.bed`
- `participants.tsv`

## Example

```bash
python scripts/1_count_snoRNA_multi_variant_carriers.py \
  --gtf annotations.gtf.gz \
  --shard-bed biallelic_shards.bed \
  --vcf-root vcfs \
  --participant-tsv participants.tsv \
  --participant-id-col platekey \
  --cpus 16 \
  --out-prefix outputs/snorna_biallelic
```

This writes:

- `outputs/snorna_biallelic.genes/0001_*.tsv`
- `outputs/snorna_biallelic.participants.tsv`
- `outputs/snorna_biallelic.gene_summary.tsv`
- `outputs/snorna_biallelic.shard_gene_map.tsv`
- `outputs/snorna_biallelic.summary.tsv`

To find participants with at least two rare variants in the same snoRNA, run:

```bash
python scripts/2_find_participants_with_2_rare_variants_same_snoRNA.py \
  --genes-dir outputs/snorna_biallelic.genes \
  --out outputs/snorna_biallelic.two_rare_same_snoRNA.tsv
```

To find participants with rare homozygous snoRNA variants, run:

```bash
python scripts/4_find_participants_with_rare_homozygous_variants.py \
  --genes-dir outputs/snorna_biallelic.genes \
  --out outputs/snorna_biallelic.homozygous_variants.tsv
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

To calculate a coverage score per snoRNA from the AGGV3 site-QC VCFs, run:

```bash
python scripts/3_calculate_snoRNA_coverage_score.py \
  --gene-summary outputs/snorna_biallelic.gene_summary.tsv \
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

If you run script 1 and your mounted directory structure differs from the default `shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz` pattern, pass `--vcf-template` with the relative path layout that matches your session.

Script 1 prints progress as it loads genes, queues each shard, and reports each shard as it finishes.
It prefers `cyvcf2` for indexed region fetches when available, then falls back to `tabix -h`, and finally to a plain scan of the VCF file. If you know the VCFs are indexed and `cyvcf2` is installed, `--region-access auto` is the fastest option.

## Participant TSV

The participant table can contain any extra columns you want. The script only needs the participant ID column.

## Notes

- The script runs without extra dependencies, but `cyvcf2` is used automatically when available.
- It treats a participant as a carrier for a site if any genotype allele is non-reference.
- Multi-allelic VCF rows are counted once per gene and participant.
- Final participant and gene summaries include only variants with `AF < 0.005`.
