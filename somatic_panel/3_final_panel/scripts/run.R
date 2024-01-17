#!/usr/bin/env Rscript
library(dplyr)
library(readr)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/somatic_panel/3_final_panel/')


# Read in tables
d.core <- readr::read_tsv('../2_core_panel/core_panel.tsv', col_types='cccccccll')
d.include <- readr::read_tsv('../resources/curation_team/genes_include.tsv', col_types='cccccccll')
d.exclude <- readr::read_tsv('../resources/curation_team/genes_exclude.tsv', col_types='cccc')


# Apply changes
d.final <- d.core |>
  dplyr::filter(! hgnc_id %in% d.exclude$hgnc_id) |>
  dplyr::bind_rows(d.include) |>
  dplyr::arrange(ensembl_gene_symbol)


# Write to disk
readr::write_tsv(d.final, 'final_panel.tsv')
