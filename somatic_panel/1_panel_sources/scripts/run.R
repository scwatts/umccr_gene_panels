#!/usr/bin/env Rscript
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(stringr)
library(tibble)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/somatic_panel/1_panel_sources/')


source('../../scripts/util.R')


# Read in HGNC data
# util.R::read_hgnc_es105
ensembl_105 <- read_ensembl_105()
# util.R::read_hgnc_latest
hgnc_latest <- read_hgnc_latest()


# Set gene list source directories and execute prep code
source_dirnames <- tibble::lst(
  dynamic=tibble::lst(
    'bushmanlab',
    'cancermine',
    'civic',
    'hmf',
    'intogen',
    'ncg',
    'oncoenrichr',
    'oncokb',
  ),
  static=tibble::lst(
    'icgc',
    'oncovar',
    'sinkala',
    'tcga',
  ),
)

# Run each prepare.R script and collect relevant data
gene_data <- c(
  purrr::map(source_dirnames$dynamic, \(x) list(n=x, p=fs::path('sources_dynamic', x))),
  purrr::map(source_dirnames$static, \(x) list(n=x, p=fs::path('sources_static', x)))
) |>
  purrr::map(\(x) {
    PREFIX <- x$p
    source(fs::path(PREFIX, 'scripts', 'prepare.R'), local=TRUE) |>
      purrr::pluck('value') |>
      dplyr::mutate(data_source=x$n) |>
      dplyr::select(dplyr::any_of(c('hgnc_id', 'data_source', 'oncogene', 'tsgene')))
  }) |>
  dplyr::bind_rows()


# Add Ensembl gene IDs and symbols
# The relationship between HGNC data and Ensembl 105 is a one-to-multiple relationship, hence can introduce additional records. Futher, since this is
# a combined list the relationship then becomes many-to-many
gene_data <- gene_data |>
  dplyr::left_join(
    dplyr::select(ensembl_105, hgnc_id, ensembl_gene_id, ensembl_transcript_id, symbol),
    by='hgnc_id',
    relationship='many-to-many',
  ) |>
  dplyr::mutate(ensembl_gene_id=stringr::str_remove(ensembl_gene_id, '\\.[0-9]+$')) |>
  dplyr::rename(ensembl_gene_symbol=symbol)


# Add in HGNC gene symbol, not all entries have an Ensembl gene ID
gene_data <- gene_data |>
  dplyr::left_join(
    dplyr::select(hgnc_latest$canonical, hgnc_id, symbol),
    by='hgnc_id',
  ) |>
  dplyr::rename(hgnc_symbol=symbol)


# Order columns and rows
gene_data <- gene_data |>
  dplyr::relocate(ensembl_gene_symbol, ensembl_gene_id, hgnc_id, hgnc_symbol, ensembl_transcript_id, oncogene, tsgene, data_source) |>
  dplyr::arrange(data_source, ensembl_gene_symbol)


# Write to disk
readr::write_tsv(gene_data, 'panel_source_data.tsv')
