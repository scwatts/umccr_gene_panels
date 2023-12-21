library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in tables
d <- dplyr::bind_rows(
  readr::read_delim(
    fs::path(PREFIX, 'data', 'ICGC.PanCancer.onco.genes.OncoVar.tsv'),
    col_select=c('Gene_symbol', 'Ensembl_ID'),
    col_types='cc',
  ),
  readr::read_delim(
    fs::path(PREFIX, 'data', 'TCGA.PanCancer.onco.genes.OncoVar.tsv'),
    col_select=c('Gene_symbol', 'Ensembl_ID'),
    col_types='cc',
  ),
) |>
  dplyr::distinct()


# Match against hgnc_es105
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_es105, by=c('Gene_symbol'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols (previous)
d.m <- d |>
  dplyr::filter(! Gene_symbol %in% d.m$Gene_symbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('Gene_symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Check provided and assigned Ensembl gene IDs
# These look okay at cursory review
# v <- apply(d.m, 1, function(r) {
#   r[['Ensembl_ID']] == stringr::str_remove(r[['data']][['ensembl_gene_id']], '\\.[0-9]+$')
# })
# View(d.m[!v, ])


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(Gene_symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
