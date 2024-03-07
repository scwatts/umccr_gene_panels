#!/usr/bin/env Rscript
library(dplyr)
library(purrr)
library(readr)
library(tibble)


#setwd('/Users/stephen/repos/gene_panels/somatic_panel/2_core_panel/')


source('../../scripts/util.R')


# Read in table
d <- readr::read_delim(
  '../1_panel_sources/panel_source_data.tsv',
  col_types='ccccccllc',
)


# Keep all OncoKB and HMF genes, split others for further processing
# Set grouping
d.grouped.retain <- d |>
  dplyr::group_by(
    retain=dplyr::if_else(data_source %in% c('oncokb', 'hmf'), 'retain', 'other')
  )

# Split groups, set list keys
v.retain_groups <- d.grouped.retain |>
  dplyr::group_split(.keep=FALSE) |>
  purrr::set_names(dplyr::group_keys(d.grouped.retain) |> dplyr::pull(1))

# Create list of HGNC IDs to retain
v.retain <- v.retain_groups$retain |>
  dplyr::pull(hgnc_id) |>
  unique()

# Set other gene entries to process
d.process <- v.retain_groups$other |>
  dplyr::filter(! hgnc_id %in% v.retain)


# Process other entries
# First require entries to be present across at least n lists
# Plot distribution to decide empirically
if (FALSE) {
  library(ggplot2)
  library(patchwork)
  library(scales)

  threshold <- 3

  # Get frequency data
  d.process.counts <- d.process |>
    dplyr::group_by(hgnc_id) |>
    dplyr::summarise(n=n())

  # Get cumulative counts
  d.c <- tibble::tibble(
    i=1:max(d.process.counts$n),
    c=sapply(1:max(d.process.counts$n), function(i) { sum(d.process.counts$n >= i) }),
  )

  # Plot frequency
  g.d <- ggplot2::ggplot(d.process.counts, ggplot2::aes(x=n))
  g.d <- g.d + ggplot2::geom_freqpoly(binwidth=1)
  g.d <- g.d + ggplot2::geom_vline(xintercept=3, colour='red')
  g.d <- g.d + ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(9), limits=c(1, 9))
  g.d <- g.d + ggplot2::scale_y_continuous(breaks=scales::pretty_breaks(5), limits=c(NA, 5000))
  g.d <- g.d + ggplot2::theme_bw()
  g.d <- g.d + ggplot2::ggtitle('Frequency distribution')

  # Plot cumulative
  g.c <- ggplot2::ggplot(d.c, ggplot2::aes(x=i, y=c))
  g.c <- g.c + ggplot2::geom_line()
  g.c <- g.c + ggplot2::geom_vline(xintercept=3, colour='red')
  g.c <- g.c + ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(9), limits=c(1, 9))
  g.c <- g.c + ggplot2::theme_bw()
  g.c <- g.c + ggplot2::ggtitle('Inverted cumulative distribution (sum records where fq ≥ i)')
  g.c <- g.c + ggplot2::labs(y='count')

  # Combine
  g.d / g.c
}


# Apply record frequency threshold; selected threshold of 3 as this follows the distribution elbow
v.retain <- d.process |>
  dplyr::group_by(hgnc_id) |>
  dplyr::summarise(n=n()) |>
  dplyr::filter(n >= 3) |>
  dplyr::pull(hgnc_id) |>
  c(v.retain)


# Create table with selected genes and set role
d.retain <- d |>
  dplyr::filter(hgnc_id %in% v.retain) |>
  dplyr::group_by(hgnc_id, ensembl_gene_id, ncbi_gene_id) |>
  dplyr::summarise(
    ensembl_gene_symbol=unique(ensembl_gene_symbol),
    hgnc_symbol=unique(hgnc_symbol),
    refseq_gene_symbol=unique(refseq_gene_symbol),
    # util.R::set_gene_role
    oncogene=set_gene_role(dplyr::across(dplyr::everything()), 'oncogene'),
    # # util.R::set_gene_role
    tsgene=set_gene_role(dplyr::across(dplyr::everything()), 'tsgene'),
  )


# Order columns and rows
d.retain <- d.retain |>
  dplyr::relocate(ensembl_gene_symbol, ensembl_gene_id, hgnc_symbol, hgnc_id, refseq_gene_symbol, ncbi_gene_id, oncogene, tsgene) |>
  dplyr::arrange(ensembl_gene_symbol)


# Write to disk
readr::write_tsv(d.retain, 'core_panel.new.tsv')
