#!/usr/bin/env Rscript
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(tibble)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/germline_panel/2_final_panel/')


# Read in table
d <- readr::read_delim(
  '../1_panel_sources/panel_source_data.tsv',
  col_types='ccccccc',
)


# Keep all genes
d.retain <- d |>
  dplyr::select(-data_source) |>
  dplyr::distinct()


# Order rows
d.retain <- d.retain |>
  dplyr::arrange(ensembl_gene_symbol)


# Write
readr::write_tsv(d.retain, 'final_panel.tsv')
