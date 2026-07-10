# snoRNA_biallelic

Minimal workflow for finding individuals with multiple snoRNA variants in AGGV3 VCF shards.

## What this repo does

- Reads a GTF and builds snoRNA intervals from `gene_type` or `gene_biotype`.
- Scans one directory of VCF or VCF.GZ shards.
- Counts snoRNA-overlapping variant sites per participant.
- Optionally joins a participant annotation TSV so case/control comparisons can happen downstream.

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
  --out-prefix outputs/snorna_biallelic
```

This writes:

- `outputs/snorna_biallelic.participants.tsv`
- `outputs/snorna_biallelic.summary.tsv`

If your mounted directory structure differs from the default `shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz` pattern, pass `--vcf-template` with the relative path layout that matches your session.

## Participant TSV

The participant table can contain any extra columns you want. The script only needs the participant ID column and, if you want grouped summaries, a case/control column.

## Notes

- The script is intentionally standard-library only.
- It treats a participant as a carrier for a site if any genotype allele is non-reference.
- Multi-allelic VCF rows are counted once per participant per site.
