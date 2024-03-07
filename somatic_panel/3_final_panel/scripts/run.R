#!/usr/bin/env Rscript
library(dplyr)
library(readr)


#setwd('/Users/stephen/repos/gene_panels/somatic_panel/3_final_panel/')


# Read in tables
d.core <- readr::read_tsv('../2_core_panel/core_panel.tsv', col_types='ccccccll')
d.include <- readr::read_tsv('../resources/curation_team/genes_include.tsv', col_types='ccccccll')
d.exclude <- readr::read_tsv('../resources/curation_team/genes_exclude.tsv', col_types='cccc')


# Apply changes
d.final <- d.core |>
  dplyr::filter(! hgnc_id %in% d.exclude$hgnc_id) |>
  dplyr::bind_rows(d.include) |>
  dplyr::arrange(ensembl_gene_symbol)


# Write to disk
readr::write_tsv(d.final, 'final_panel.tsv')
