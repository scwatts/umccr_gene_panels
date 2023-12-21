#!/usr/bin/env Rscript
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(stringr)
library(tibble)


#setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/somatic_panel/1_sources/')


source('../../scripts/util.R')


# Read in HGNC data
# util.R::read_hgnc_es105
hgnc_es105 <- read_hgnc_es105()
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


# Add HGNC data (HGNC ID, Ensembl gene ID)
# Match against hgnc_es105
gene_data.matched <- gene_data |>
  dplyr::nest_join(hgnc_es105, by='hgnc_id', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  tidyr::hoist(data, 'ensembl_gene_id') |>
  dplyr::mutate(hgnc_source='es105')

# Match against hgnc_latest for missing entries
# Some genes do not have an Ensembl gene ID
gene_data.matched <- gene_data |>
  dplyr::filter(! hgnc_id %in% gene_data.matched$hgnc_id) |>
  dplyr::nest_join(hgnc_latest$canonical, by='hgnc_id', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  tidyr::hoist(data, 'ensembl_gene_id') |>
  dplyr::mutate(hgnc_source='latest') |>
  dplyr::bind_rows(gene_data.matched)


# Add HGNC gene symbol, sort rows, drop and order cols
gene_data <- gene_data.matched |>
  dplyr::mutate(ensembl_gene_id=stringr::str_remove(ensembl_gene_id, '\\.[0-9]+$')) |>
  tidyr::hoist(data, 'symbol') |>
  dplyr::arrange(data_source, symbol) |>
  dplyr::select(-data) |>
  dplyr::relocate(
    hgnc_id,
    symbol,
    ensembl_gene_id,
    hgnc_source,
    oncogene,
    tsgene,
    data_source,
  )


# Write to disk
readr::write_tsv(gene_data, 'prepared_combined.tsv')
