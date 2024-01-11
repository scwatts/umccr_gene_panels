library(dplyr)
library(fs)
library(readr)
library(stringr)
library(tidyr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'NCG_cancerdrivers_annotation_supporting_evidence.tsv'),
  col_select=c('symbol', 'pubmed_id', 'cgc_annotation', 'vogelstein_annotation', 'saito_annotation', 'NCG_oncogene', 'NCG_tsg'),
  col_types='ccccc',
)


# Determine gene role, allowing any annotation to set status affirmatively
d.p <- d |>
  dplyr::group_by(symbol) |>
  dplyr::summarize(
    oncogene=(
      stringr::str_equal('1', NCG_oncogene) |
        stringr::str_detect(cgc_annotation, 'oncogene') |
        stringr::str_detect(vogelstein_annotation, 'Oncogene') |
        stringr::str_detect(saito_annotation, 'Oncogene')
    ) |> any() |> tidyr::replace_na(FALSE),
    tsgene=(
      stringr::str_equal('1', NCG_tsg) |
        stringr::str_detect(cgc_annotation, 'TSG') |
        stringr::str_detect(vogelstein_annotation, 'TSG') |
        stringr::str_detect(saito_annotation, 'TSG')
    ) |> any() |> tidyr::replace_na(FALSE),
  )


# Match against HGNC latest
# Symbols
d.m <- d.p |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))

# Symbols (previous)
d.m <- d.p |>
  dplyr::filter(! symbol %in% d.m$symbol) |>
  dplyr::nest_join(hgnc_latest$previous, by=c('symbol'='prev_symbol'), name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1)) |>
  dplyr::bind_rows(d.m)


# Hoist matched HGNC from nested tibble, then format
d.m <- d.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)


# Write to disk
readr::write_tsv(d.m, file=fs::path(PREFIX, 'prepared.tsv'))
