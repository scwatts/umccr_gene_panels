library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'gene-ids-for-set-Cancer Gene Census.tsv'),
  col_names=c('ensembl_gene_id', 'symbol'),
  col_types='cc',
)


# Match against hgnc_es105
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_es105, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Symbols (previous); completes all recording matches
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Check provided and assigned Ensembl gene IDs
# Each non-match manually checked via ICGC and genenames.org:
#   - DUX4, HIST1H3B, HIST1H4I, MLLT6, NCOA4, SDHD, SSX4, TAF15
# v <- apply(d.m, 1, function(r) {
#   r[['ensembl_gene_id']] == stringr::str_remove(r[['data']][['ensembl_gene_id']], '\\.[0-9]+$')
# })
# d.m[!v, ]


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
