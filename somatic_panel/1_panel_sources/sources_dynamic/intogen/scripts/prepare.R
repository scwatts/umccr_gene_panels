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


# Match against HGNC latest
# Symbols
d.m <- d.u |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Hoist matched HGNC from nested tibble
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
