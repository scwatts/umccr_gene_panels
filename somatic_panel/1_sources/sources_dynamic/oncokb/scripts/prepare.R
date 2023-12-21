library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'cancerGeneList.txt'),
  col_select=c(
    'symbol'='Hugo Symbol',
    'oncogene'='Is Oncogene',
    'tsgene'='Is Tumor Suppressor Gene',
    'oncokb_annotated'='OncoKB Annotated'
  ),
  col_types='cccc',
)


# Select only OncoKB annotated entries
d.f <- d |>
  dplyr::filter(oncokb_annotated == 'Yes')


# Match against hgnc_es105
# Symbols
d.m <- d.f |>
  dplyr::nest_join(hgnc_es105, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols
# H3P6: processed_pseudogene; excluded in processing of es105
# TRA, TRD, TRD, TRG: not present in es105
d.m <- d.f |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Set oncogene and tsgene, remove other unneeded column
d.p <- d.m |>
  dplyr::mutate(
    oncogene=oncogene == 'Yes',
    tsgene=tsgene == 'Yes',
  ) |>
  dplyr::select(-oncokb_annotated)


# Hoist matched HGNC from nested tibble
d.p <- d.p |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.p, file=fs::path(PREFIX, 'prepared.tsv'))
