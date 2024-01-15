library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(fs::path(PREFIX, 'data', '01-Dec-2023-GeneSummaries.tsv'), col_types='ccccccc')


# Apply fixes for matching
# COX2 matches on ambiguous alias symbol in hgnc_latest; manually confirmed to be MT-CO2
d$name[d$name=='COX2'] <- 'MT-CO2'
# ND1 matches on ambiguous alias symbol in hgnc_latest; manually confirmed to be MT-ND1
d$name[d$name=='ND1'] <- 'MT-ND1'


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('name'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(name)



# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))


# Format consistent processing in sourcing script
d.m <- d.m |>
  dplyr::rename(symbol=name)
