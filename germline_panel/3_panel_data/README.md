# Panel data

Generate panel data for UMCCR post-processing

```bash
mkdir -p output/

cut -f2 -d$'\t' ../2_final_panel/final_panel.tsv | \
  tail -n+2 | \
  grep -v NA | \
  sed 's/\..*$//' | \
  sort > output/umccr_predispoition_genes.latest.genes

./scripts/create_subset_bed.py > output/umccr_predispoition_genes.cds.bed \
  --panel_fp ../2_final_panel/final_panel.tsv \
  --ensembl_cds_data_fp ../../resources/ensembl_cds_data.tsv
```
