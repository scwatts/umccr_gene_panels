library(assertthat)
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(tidyr)


# Read in table and ensure no duplicate symbols
d <- readr::read_delim(fs::path(PREFIX, 'data', 'allOnco_June2021.tsv'), col_types='ccccc')
assertthat::assert_that(all(table(d$symbol) == 1))


# Apply fixes for matching
# Correct letter casing
d$symbol[d$symbol=='C15ORF65'] <- 'C15orf65'
# GGTA1P has two matches, one being the inactive gene; rename to ensure we match the correct one only
d$symbol[d$symbol=='GGTA1P'] <- 'GGTA2P'


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Symbols (previous)
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  # NOTE(SW): GGTA1P has two Ensembl gene entries, must unnest
  tidyr::unnest_longer(hgnc_id) |>
  dplyr::select(-data) |>
  dplyr::arrange(as.numeric(row))


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
