#!/usr/bin/env Rscript


library(dplyr)
library(purrr)
library(readr)
library(tidyr)


# setwd('/Users/stephen/projects/umccr_key_genes_202312/2_sources/repo/fusion_panel/1_panel_data/')


source('../../scripts/util.R')


# Read in HGNC data
# util.R::read_hgnc_es105
ensembl_105 <- read_ensembl_105()
# util.R::read_hgnc_latest
hgnc_latest <- read_hgnc_latest()


# Read in data
d.hmf <- readr::read_delim('data/known_fusion_data.38.csv', col_types='c')
d.somatic_panel <- readr::read_delim('../../somatic_panel/3_final_panel/final_panel.tsv', col_types='ccccll')


# Set gene symbol for somatic panel
# NOTE(SW): some entries only have HGNC symbol, so use that if necessary
d.somatic_panel <- d.somatic_panel |>
  dplyr::mutate(
    symbol=dplyr::case_when(
      !is.na(ensembl_gene_symbol) ~ ensembl_gene_symbol,
      !is.na(hgnc_symbol) ~ hgnc_symbol,
    )
  )


# Get dataframe of gene symbols for matching
d.hmf.genes <- tibble::tibble(
  symbol=d.hmf |>
    dplyr::select(FiveGene, ThreeGene) |>
    unlist() |>
    unique() |>
    purrr::discard(is.na)
)


# Match against HGNC latest
# Symbols
d.hmf.genes.m <- d.hmf.genes |>
  dplyr::nest_join(hgnc_latest$canonical, by='symbol', name='data', keep=TRUE) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) >= 1))


# NOTE(SW): AL121790.1 and C19MC are allowed to be dropped here


# Hoist matched HGNC from nested tibble, then format
d.hmf.genes <- d.hmf.genes.m |>
  tidyr::hoist(data, 'hgnc_id') |>
  dplyr::select(-data) |>
  dplyr::arrange(symbol)

d.hmf <- d.hmf |>
  dplyr::left_join(
    dplyr::rename(d.hmf.genes, FiveGene_hgnc_id=hgnc_id),
    by=c('FiveGene'='symbol'),
  ) |>
  dplyr::left_join(
    dplyr::rename(d.hmf.genes, ThreeGene_hgnc_id=hgnc_id),
    by=c('ThreeGene'='symbol'),
  )


# Remove entries that are not part of the somatic gene panel
# This currently only excludes C19MC, which was dropped above since it has no HGNC or Ensembl entry
d.hmf <- d.hmf |>
  dplyr::filter(
    FiveGene_hgnc_id %in% d.hmf.genes$hgnc_id | ThreeGene_hgnc_id %in% d.hmf.genes$hgnc_id,
  )


# Set aside dataframe to write out
d.hmf.subset <- d.hmf |>
  dplyr::select(-c(FiveGene_hgnc_id, ThreeGene_hgnc_id))


# Add entries from the somatic panel to the fusion panel with the promiscuous fusion type
# First, get lists of genes that do not have a promiscuous entry already
# NOTE(SW): still allowing PROMISCUOUS_5 for IGH, IGK despite these already having IG_PROMISCUOUS since that is restricted
v.promiscuous_types <- tibble::lst(
  five='PROMISCUOUS_5',
  three='PROMISCUOUS_3',
)

v.genes_add <- lapply(v.promiscuous_types, \(x) {
  d.hmf |>
    dplyr::filter(Type==x) |>
    dplyr::pull(
      dplyr::case_when(
        x == 'PROMISCUOUS_5' ~ FiveGene_hgnc_id,
        x == 'PROMISCUOUS_3' ~ ThreeGene_hgnc_id,
      )
    ) |>
    unlist() |>
    unique() |>
    purrr::discard(is.na)
})

# Create partial rows for new entries
d.hmf.new_data.partial <- dplyr::bind_rows(
  tibble::tibble(
    Type='PROMISCUOUS_5',
    FiveGene=d.somatic_panel |> dplyr::filter(! hgnc_id %in% v.genes_add$five) |> dplyr::pull(symbol),
    ThreeGene=NA,
  ),
  tibble::tibble(
    Type='PROMISCUOUS_3',
    FiveGene=NA,
    ThreeGene=d.somatic_panel |> dplyr::filter(! hgnc_id %in% v.genes_add$three) |> dplyr::pull(symbol),
  ),
)

# Complete rows of each entry
d.hmf.new_data <- d.hmf.new_data.partial |>
  dplyr::bind_cols(
    tibble::tibble(
      CancerTypes=NA,
      PubMedId=NA,
      KnownExonTranscript=NA,
      KnownExonUpRange=NA,
      KnownExonDownRange=NA,
      HighImpactPromiscuous=NA,
      Overrides=NA,
    )
  ) |>
  # Sort more nicely (keeping new five/three promiscuous entries together)
  dplyr::mutate(
    sort=dplyr::case_when(
      !is.na(FiveGene) ~ FiveGene,
      !is.na(ThreeGene) ~ ThreeGene,
    ),
  ) |>
  dplyr::arrange(sort) |>
  dplyr::select(-sort)


# Combine modified Hartwig fusion data with new entries
d.hmf <- d.hmf |>
  dplyr::bind_rows(d.hmf.new_data) |>
  dplyr::select(-c(FiveGene_hgnc_id, ThreeGene_hgnc_id))


# Finally add two extra columns + rename another to make compatible with gene-utils from hmftools
d.hmf.database <- d.hmf |>
  dplyr::mutate(
    OverridesRef38=Overrides,
    Overrides=NA,
    KnownExonTranscriptRef38=NA,
  )


# Write to disk
readr::write_csv(d.hmf.subset, file='output/hmftools_fusion_data.subset.csv', na='')
readr::write_tsv(d.hmf.database, file='output/fusion_database.tsv', na='')
