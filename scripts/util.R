library(dplyr)
library(here)
library(readr)
library(stringr)
library(tidyr)


read_ensembl_105 <- function() {
  fp <- here::here('resources/ensembl.genes.tsv')
  readr::read_delim(fp, delim='\t', col_types='c') |>
    # Exclude genes with no HGNC ID
    dplyr::filter(!is.na(hgnc_id)) |>
    # util.R::exclude_genes_in_par
    exclude_genes_in_par() |>
    # Remove ENSG00000269900.3 since it is a duplicate gene in 105 that is removed in the latest
    dplyr::filter(gene_id != 'ENSG00000269900.3') |>
    # Make gene ID column name more specific
    dplyr::rename(ensembl_gene_id=gene_id)
}


read_refseq <- function() {
  fp <- here::here('resources/refseq.genes.tsv')
  readr::read_delim(fp, delim='\t', col_types='c') |>
    # Exclude genes with no HGNC ID
    dplyr::filter(!is.na(hgnc_id)) |>
    # util.R::exclude_genes_in_par
    exclude_genes_in_par() |>
    # Make gene ID column name more specific
    dplyr::rename(ncbi_gene_id=gene_id)
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


exclude_genes_in_par <- function(.data) {
  # Regions from https://www.ensembl.org/info/genome/genebuild/human_PARS.html

  # NOTE(SW): this additionally filters XGY2 (HGNC:34022), which spans the PAR boundary but is for some reason not considered a PAR gene; ignoring since it is not an relevant gene

  .data |>
    # PAR1 on chrY
    dplyr::filter(
      ! (contig == 'chrY' & (end >= 10001 & start <= 2781479))
    ) |>
    # PAR2 on chrY
    dplyr::filter(
      ! (contig == 'chrY' & (end >= 56887903 & start <= 57217415))
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
