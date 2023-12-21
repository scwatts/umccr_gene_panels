library(dplyr)
library(fs)
library(readr)
library(stringr)
library(tidyr)


# Read in table
d <- readr::read_csv(
  fs::path(PREFIX, 'data', 'table_s1.csv'),
  col_types='cccccccccccc',
)


# Set gene role
d.p <- d |>
  dplyr::rename(gene_role=`Tumor suppressor or oncogene prediction (by 20/20+)`) |>
  dplyr::group_by(Gene) |>
  dplyr::summarize(
    oncogene=stringr::str_detect(gene_role, 'oncogene') |> any() |> tidyr::replace_na(FALSE),
    tsgene=stringr::str_detect(gene_role, 'tsg') |> any() |> tidyr::replace_na(FALSE),
  )


# Match against hgnc_es105
# Symbols
d.m <- d.p |>
  dplyr::nest_join(hgnc_es105, by=c('Gene'='symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# Match against hgnc_latest
# Symbols (previous)
d.m <- d.p |>
  dplyr::filter(! Gene %in% d.m$Gene) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('Gene'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(Gene)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
