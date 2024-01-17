library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_csv(
  fs::path(PREFIX, 'data', 'cancer_genes.csv'),
  col_types='cccc',
)


# Apply fixes for matching
# ZEB1 has two entries, selecting the first
d <- d |>
  dplyr::filter(
    dplyr::row_number() != which(d$HugoSymbol == 'ZEB1') |> tail(1)
  )


# Set gene role
d.p <- d |>
  dplyr::mutate(
    oncogene=RoleInCancer == 'Oncogenes',
    tsgene=RoleInCancer == 'TSGs',
  )


# Match against HGNC latest
# Symbols
d.m <- d.p |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('HugoSymbol'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Symbols (previous)
d.m <- d.p |>
  dplyr::filter(! HugoSymbol %in% d.m$HugoSymbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('HugoSymbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(HugoSymbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))


# Format consistent processing in sourcing script
d.m <- d.m |>
  dplyr::rename(symbol=HugoSymbol)
