#!/usr/bin/env Rscript
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(tibble)


#setwd('/Users/stephen/repos/gene_panels/somatic_panel/1_panel_sources/')


source('../../scripts/util.R')


# Read in HGNC data
# util.R::read_hgnc_es105
ensembl_105 <- read_ensembl_105()
# util.R::read_refseq
refseq <- read_refseq()
# util.R::read_hgnc_latest
hgnc_latest <- read_hgnc_latest()


# Set gene list source directories and execute prep code
source_dirnames <- tibble::lst(
  dynamic=tibble::lst(
    'bushmanlab',
    'cancermine',
    'civic',
    'cpsr',
    'hmf',
    'intogen',
    'ncg',
    'oncoenrichr',
    'oncokb',
    'pmcc_fcc',
  ),
  static=tibble::lst(
    'icgc',
    'oncovar',
    'sinkala',
    'tcga',
    'tso500',
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


# Add Ensembl 105 gene IDs and symbols, not all entries have Ensembl records
gene_data <- gene_data |>
  dplyr::left_join(
    dplyr::select(ensembl_105, hgnc_id, ensembl_gene_id, symbol),
    by='hgnc_id',
  ) |>
  dplyr::rename(ensembl_gene_symbol=symbol)


# Add in latest HGNC gene symbol
gene_data <- gene_data |>
  dplyr::left_join(
    dplyr::select(hgnc_latest$canonical, hgnc_id, symbol),
    by='hgnc_id',
  ) |>
  dplyr::rename(hgnc_symbol=symbol)


# Add latest RefSeq/NCBI data
gene_data <- gene_data |>
  dplyr::left_join(
    dplyr::select(refseq, hgnc_id, ncbi_gene_id, symbol),
    by='hgnc_id',
  ) |>
  dplyr::rename(refseq_gene_symbol=symbol)

# NOTE(SW): HGNC:4641 (GSTT1) is currently on chr22 ALT in GRCh38, hence only has HGNC identifiers

# Order columns and rows
gene_data <- gene_data |>
  dplyr::relocate(ensembl_gene_symbol, ensembl_gene_id, hgnc_symbol, hgnc_id, refseq_gene_symbol, ncbi_gene_id, oncogene, tsgene, data_source) |>
  dplyr::arrange(data_source, ensembl_gene_symbol)


# Write to disk
readr::write_tsv(gene_data, 'panel_source_data.tsv')
