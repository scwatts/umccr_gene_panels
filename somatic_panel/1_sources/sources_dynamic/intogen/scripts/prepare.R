library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', '2023-05-31_IntOGen-Drivers', 'Compendium_Cancer_Genes.tsv'),
  col_select=c('symbol'='SYMBOL'),
  col_types='c',
)


# Remove duplicate symbol entries
d.u <- d |>
  dplyr::distinct()


# Match against hgnc_es105
# Symbols
d.m <- d.u |>
  dplyr::nest_join(hgnc_es105, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols (the single entry here, SLC35E2A, is in es105 but excluded during processing since it is set as `transcribed_unprocessed_pseudogene`)
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
