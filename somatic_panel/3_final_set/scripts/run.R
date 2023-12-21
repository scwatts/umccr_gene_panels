#!/usr/bin/env Rscript
library(dplyr)
library(readr)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/somatic_panel/3_final_set/')


# Read in tables
d.core <- readr::read_tsv('../2_full_set/gene_list.tsv', col_types='cccll')
d.include <- readr::read_tsv('../resources/curation_team/genes_include.tsv', col_types='cccll')
d.exclude <- readr::read_tsv('../resources/curation_team/genes_exclude.tsv', col_types='ccc')


# Apply changes
d.final <- d.core |>
  dplyr::filter(! hgnc_id %in% d.exclude$hgnc_id) |>
  dplyr::bind_rows(d.include)


# Write to disk
readr::write_tsv(d.final, 'gene_list.tsv')
