library(assertthat)
library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'cpsr_superpanel_2022_01.tsv'),
  col_select=c('symbol', 'ensembl_gene_id'),
  col_types='cc',
)


# Apply fixes
# Update RMRP (HGNC:10031) Ensembl gene ID to expect record to pass checks below
d$ensembl_gene_id[d$ensembl_gene_id=='ENSG00000269900'] <- 'ENSG00000277027'
# Set GBA to canonical GBA1 symbol to simplify matching below
d$symbol[d$symbol=='GBA'] <- 'GBA1'


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id', 'hgnc_ensembl_gene_id'='ensembl_gene_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Check we've matched correctly
assertthat::assert_that(all(d.m$ensembl_gene_id==d.m$hgnc_ensembl_gene_id))


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
