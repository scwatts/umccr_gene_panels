library(dplyr)
library(fs)
library(readr)
library(tibble)
library(tidyr)


# Read in table
d <- readr::read_lines(
  fs::path(PREFIX, 'data', 'pmcc_fcc.txt'),
) |>
  tibble::enframe(name=NULL, value='symbol')


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
