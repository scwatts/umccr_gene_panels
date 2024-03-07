library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in tables
# Somatic panel
d.panel <- readr::read_delim(
  fs::path(PREFIX, 'data', 'DriverGenePanel.38.tsv'),
  col_select=c('gene', 'likelihoodType'),
  col_types='cc',
)

# Fusion panel
d.fusions <- readr::read_delim(
  fs::path(PREFIX, 'data', 'known_fusion_data.38.csv'),
  col_select=c('FiveGene', 'ThreeGene'),
  col_types='cc',
) |>
  tidyr::pivot_longer(
    cols=tidyselect::everything(),
    names_to='source',
    values_to='gene',
  ) |>
  dplyr::select(-source) |>
  dplyr::filter(!is.na(gene)) |>
  dplyr::distinct() |>
  dplyr::filter(! gene %in% d.panel$gene)

# Create input data tibble
d <- dplyr::bind_rows(
  d.panel,
  d.fusions,
)


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('gene'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


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
