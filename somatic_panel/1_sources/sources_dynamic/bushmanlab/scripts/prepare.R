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
# MKL1 doesn't match hgnc_es105 on symbol and instead incorrectly matches on synonym with MAL, set correct name here
d$symbol[d$symbol=='MKL1'] <- 'MRTFA'


# Match against hgnc_es105
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_es105, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Names/descriptions
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_es105, by='name', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)

# Previous symbols
# Matches look okay by gene name
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  tidyr::separate_longer_delim(prevSymbols, delim=',') |>
  dplyr::nest_join(hgnc_es105, by=c('prevSymbols'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)

# Synonyms
# Matches manually confirmed via HGNC database lookup
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  tidyr::separate_longer_delim(synonyms, delim=',') |>
  dplyr::nest_join(hgnc_es105, by=c('synonyms'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Match remaining against hgnc_latest
# Symbols
d.m <- d |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)

# Symbols (previous); completes all recording matches
# Matches manually confirmed via HGNC database lookup
d.m <- d |>
  dplyr::filter(! symbol %in% c(d.m$symbol, d.m$symbol)) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(as.numeric(row))


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
