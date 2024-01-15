library(dplyr)
library(fs)
library(readr)
library(tidyr)


# Read in tables
d <- dplyr::bind_rows(
  readr::read_delim(
    fs::path(PREFIX, 'data', 'ICGC.PanCancer.onco.genes.OncoVar.tsv'),
    col_select=c('Gene_symbol', 'Ensembl_ID'),
    col_types='cc',
  ),
  readr::read_delim(
    fs::path(PREFIX, 'data', 'TCGA.PanCancer.onco.genes.OncoVar.tsv'),
    col_select=c('Gene_symbol', 'Ensembl_ID'),
    col_types='cc',
  ),
) |>
  dplyr::distinct()


# Match against HGNC latest
# Symbols
d.m <- d |>
  dplyr::nest_join(hgnc_latest$canonical, by=c('Gene_symbol'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Symbols (previous)
d.m <- d |>
  dplyr::filter(! Gene_symbol %in% d.m$Gene_symbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('Gene_symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(Gene_symbol)


# # Include Ensembl gene IDs where available; the relationship between HGNC data and Ensembl 105 is a one-to-multiple relationship, hence can introduce additional records
# d.e <- d.m |>
#   dplyr::left_join(dplyr::select(ensembl_105, hgnc_id, ensembl_gene_id), by='hgnc_id') |>
#   dplyr::mutate(ensembl_gene_id=stringr::str_remove(ensembl_gene_id, '\\..*$'))
#
# rr <- d.e[d.e$Ensembl_ID!=d.e$ensembl_gene_id, ]
# message(paste0('https://asia.ensembl.org/Homo_sapiens/Gene/Summary?g=', rr$Ensembl_ID, '\n'))
#
# All 36 correctly updated
# * on scaffold contig: DAXX, HLA-A, STK19, TRIM27
# * SSX1: Ensembl gene ID was for SSX2
# * SSX2: Ensembl gene ID was for SSX2B
# * SSX4: Ensembl gene ID was for SSX4B
# * all others had no entry in Ensembl 105 or HGNC


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
