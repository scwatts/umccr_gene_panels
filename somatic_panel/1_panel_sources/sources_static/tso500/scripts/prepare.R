library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_csv(
  fs::path(PREFIX, 'data', 'tso500.tsv'),
  col_types='c',
  col_names='gene',
)


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('gene'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Symbols (previous)
d.m <- d |>
  dplyr::filter(! gene %in% d.m$gene) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('gene'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(gene)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))