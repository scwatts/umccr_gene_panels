library(dplyr)
library(fs)
library(readr)
library(tidyr)
library(stringr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'gene-ids-for-set-Cancer Gene Census.tsv'),
  col_names=c('ensembl_gene_id', 'symbol'),
  col_types='cc',
)


# Apply fixes
# Remove ENSG00000255292 since it doesn't match SDHD or any gene in Ensembl 105 or latest HGNC
d <- d[d$ensembl_gene_id!='ENSG00000255292', ]


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


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# # Include Ensembl gene IDs where available; the relationship between HGNC data and Ensembl 105 is a one-to-multiple relationship, hence can introduce additional records
# d.e <- d.m |>
#   dplyr::left_join(dplyr::select(ensembl_105, hgnc_id, ensembl_gene_id), by='hgnc_id') |>
#   dplyr::mutate(ensembl_gene_id.y=stringr::str_remove(ensembl_gene_id.y, '\\..*$'))
#
# d.e[d.e$ensembl_gene_id.x!=d.e$ensembl_gene_id.y, ]

# Manually confirmed incoming data is correct for the following mismatched Ensembl gene IDs; the source Ensembl gene IDs were not present in Ensembl 105 or latest HGNC
# # A tibble: 7 × 4
# ensembl_gene_id.x symbol   hgnc_id    ensembl_gene_id.y
# <chr>             <chr>    <chr>      <chr>
# 1 ENSG00000258389   DUX4     HGNC:50800 ENSG00000260596
# 2 ENSG00000124693   HIST1H3B HGNC:4776  ENSG00000286522
# 3 ENSG00000198339   HIST1H4I HGNC:4793  ENSG00000276180
# 4 ENSG00000108292   MLLT6    HGNC:7138  ENSG00000275023
# 5 ENSG00000138293   NCOA4    HGNC:7671  ENSG00000266412
# 6 ENSG00000204645   SSX4     HGNC:11338 ENSG00000268009
# 7 ENSG00000172660   TAF15    HGNC:11547 ENSG00000270647


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
