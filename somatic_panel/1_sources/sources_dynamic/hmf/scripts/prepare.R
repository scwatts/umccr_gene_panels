library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'DriverGenePanel.38.tsv'),
  col_select=c('gene', 'likelihoodType'),
  col_types='cc',
)


# Match against hgnc_es105
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_es105, by=c('gene'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols
d.m <- d |>
  dplyr::filter(! gene %in% d.m$gene) |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('gene'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Set gene role
d.m <- d.m |>
  dplyr::mutate(
    oncogene=likelihoodType == 'ONCO',
    tsgene=likelihoodType == 'TSG',
  )


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(gene)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
