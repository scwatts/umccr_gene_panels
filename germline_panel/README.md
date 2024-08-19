# Germline panel

The UMCCR germline panel is currently composed only of entries from the PMCC FCC panel with no modifications applied by
the curation team. While there is only one source list, the panel creation process is structured similarly to the
somatic panel so that additional sources may be added in the future.

## Procedure

1. Review literature, community for new sources
  * Include new sources in `1_panel_sources/` if appropriate with feedback from curation team
2. Apply any upstream changes for sources in `1_panel_sources/`
3. Construct final panel in `2_final_panel/`
4. Generate panel data files in `3_panel_data/`

## CPSR note

The following genes are absent from the CPSR superpanel (panel 0) [v20220203;
`data/grch38/virtual_panels/cpsr_superpanel.grch38.tsv`] and hence do not appear in the CPSR report:

```text
ensembl_gene_symbol  ensembl_gene_id     hgnc_symbol  hgnc_id     refseq_gene_symbol  ncbi_gene_id
CSDE1                ENSG00000009307.17  CSDE1        HGNC:29905  CSDE1               7812
EGLN1                ENSG00000135766.9   EGLN1        HGNC:1232   EGLN1               54583
EGLN2                ENSG00000269858.6   EGLN2        HGNC:14660  EGLN2               112398
```

While variants in these genes do not appear in the CPSR report, they are present in the output germline small variant
VCF file.
