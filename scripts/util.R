library(dplyr)
library(here)
library(readr)
library(tidyr)


read_hgnc_es105 <- function() {
  fp <- here::here('resources/ensembl_105_genes_with_hgnc.tsv')
  readr::read_delim(fp, delim='\t', col_types='cccc')
}


read_hgnc_latest <- function() {
  fp <- here::here('resources/hgnc_complete_set_2023-11-01.tsv')
  d <- readr::read_delim(
    fp,
    delim='\t',
    col_select=c(hgnc_id, ensembl_gene_id, symbol, alias_symbol, prev_symbol, name),
    col_types='cccccc',
  )
  
  d.alias <- d |>
    tidyr::separate_longer_delim(
      alias_symbol,
      delim='|',
    )
  
  d.previous <- d |>
    tidyr::separate_longer_delim(
      prev_symbol,
      delim='|',
    )

  list(
    canonical=d,
    alias=d.alias,
    previous=d.previous
  )
}


set_gene_role <- function(.data, .col_name) {
  
  # NOTE(SW): tibble hits a critical bug on dplyr::filter, casting to data.frame fixes
  .data <- data.frame(.data)
  
  # Return NA if there are no annotations for this gene
  if (dplyr::pull(.data, .col_name) |> is.na() |> all()) {
    return(NA)
  }
  
  # Prioritise OncoKB and oncoEnrichR annotations
  if ('oncokb' %in% .data$data_source) {
    dplyr::filter(.data, data_source=='oncokb') |> dplyr::pull(.col_name)
  } else if ('oncoenrichr' %in% .data$data_source) {
    dplyr::filter(.data, data_source=='oncoenrichr') |> dplyr::pull(.col_name)
  } else {
    dplyr::pull(.data, .col_name) |> tidyr::replace_na(FALSE) |> any()
  }
}