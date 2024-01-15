#!/usr/bin/env Rscript
library(dplyr)
library(readr)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/somatic_panel/3_final_panel/')


# Read in tables
d.core <- readr::read_tsv('../2_core_panel/core_panel.tsv', col_types='cccccll')
d.include <- readr::read_tsv('../resources/curation_team/genes_include.tsv', col_types='ccccll')
d.exclude <- readr::read_tsv('../resources/curation_team/genes_exclude.tsv', col_types='cccc')


# Apply changes
d.final <- d.core |>
  dplyr::filter(! hgnc_id %in% d.exclude$hgnc_id) |>
  dplyr::bind_rows(d.include) |>
  dplyr::arrange(ensembl_gene_symbol)


# Exclude PAR region genes where there is another non-PAR entry
v.exclude_par <- d.final |>
  dplyr::filter(
    duplicated(ensembl_gene_symbol) &
      ! is.na(ensembl_gene_symbol) &
      stringr::str_detect(ensembl_gene_id, '_PAR_Y$')
  ) |>
  dplyr::pull(ensembl_gene_id)

d.final <- d.final |>
  dplyr::filter(! ensembl_gene_id %in% v.exclude_par)


# Write to disk
readr::write_tsv(d.final, 'final_panel.tsv')
