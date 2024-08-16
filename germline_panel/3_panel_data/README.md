# Panel data

Generate panel data for UMCCR post-processing

```bash
mkdir -p output/

cut -f2 -d$'\t' ../2_final_panel/final_panel.tsv | \
  tail -n+2 | \
  grep -v NA | \
  sed 's/\..*$//' | \
  sort > output/umccr_predisposition_genes.gene_symbols.txt

./scripts/create_transcript_bed.py > output/umccr_predisposition_genes.transcript_regions.bed \
  --panel_fp ../2_final_panel/final_panel.tsv \
  --ensembl_transcript_data_fp ../../resources/ensembl.transcripts.bed
```
