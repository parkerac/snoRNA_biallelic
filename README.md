# snoRNA_biallelic

Minimal workflow for finding individuals with multiple snoRNA variants in AGGV3 VCF shards.

## What this repo does

- Reads a GTF and builds snoRNA intervals from `gene_type` or `gene_biotype`.
- Uses `biallelic_shards.bed` to find the relevant shard/subshard VCFs for each snoRNA gene.
- Writes one TSV per gene as soon as that gene finishes.
- Produces end-of-run gene and participant summaries.
- Records the participant genotype in each per-gene TSV.
- Merges nearby genes into shared fetch windows so each shard is queried fewer times.

## Suggested layout on CloudOS

- `annotations.gtf.gz`
- `vcfs/`
- `biallelic_shards.bed`
- `participants.tsv`

## Example

```bash
python scripts/count_snoRNA_multi_variant_carriers.py \
  --gtf annotations.gtf.gz \
  --shard-bed biallelic_shards.bed \
  --vcf-root vcfs \
  --participant-tsv participants.tsv \
  --participant-id-col platekey \
  --group-col phenotype \
  --cpus 16 \
  --out-prefix outputs/snorna_biallelic
```

This writes:

- `outputs/snorna_biallelic.genes/0001_*.tsv`
- `outputs/snorna_biallelic.participants.tsv`
- `outputs/snorna_biallelic.gene_summary.tsv`
- `outputs/snorna_biallelic.shard_gene_map.tsv`
- `outputs/snorna_biallelic.summary.tsv`

If your mounted directory structure differs from the default `shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz` pattern, pass `--vcf-template` with the relative path layout that matches your session.

The script prints progress as it loads genes, queues each shard, and reports each shard as it finishes.
It prefers `cyvcf2` for indexed region fetches when available, then falls back to `tabix -h`, and finally to a plain scan of the VCF file. If you know the VCFs are indexed and `cyvcf2` is installed, `--region-access auto` is the fastest option.

## Participant TSV

The participant table can contain any extra columns you want. The script only needs the participant ID column and, if you want grouped summaries, a case/control column.

## Notes

- The script runs without extra dependencies, but `cyvcf2` is used automatically when available.
- It treats a participant as a carrier for a site if any genotype allele is non-reference.
- Multi-allelic VCF rows are counted once per gene and participant.
