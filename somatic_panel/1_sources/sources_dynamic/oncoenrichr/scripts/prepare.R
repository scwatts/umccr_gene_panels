library(dplyr)
library(fs)
library(readr)
library(stringr)


# Read in table
d <- readr::read_delim(
  fs::path(PREFIX, 'data', 'genedb_v1.4.2.tsv'),
  col_select=c('hgnc_id', 'tumor_suppressor', 'oncogene', 'cancer_max_rank'),
  col_type='cccd',
)


# Review cancer_max_rank to empirically decide threshold; use arbitrary 0.5 threshold, rank is fairly well scaled b/n 0-1
if (FALSE) {
  library(ggplot2)
  library(patchwork)

  threshold <- 0.5

  g.fq <- ggplot2::ggplot(d[d$cancer_max_rank>0, ], aes(x=cancer_max_rank))
  g.fq <- g.fq + ggplot2::geom_freqpoly(binwidth=0.01)
  g.fq <- g.fq + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.fq <- g.fq + ggplot2::theme_bw()
  g.fq <- g.fq + ggplot2::theme(axis.title.x=ggplot2::element_blank())
  g.fq <- g.fq + ggplot2::ggtitle('Frequency distribution')

  g.cl <- ggplot2::ggplot(d[d$cancer_max_rank>0, ], aes(x=cancer_max_rank))
  g.cl <- g.cl + ggplot2::geom_abline(slope=1, intercept=0, color='grey')
  g.cl <- g.cl + ggplot2::stat_ecdf(geom='step')
  g.cl <- g.cl + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.cl <- g.cl + ggplot2::theme_bw()
  g.cl <- g.cl + ggplot2::ggtitle('Cumulative distribution')
  g.cl <- g.cl + ggplot2::labs(x='Cancer rank (max)')

  g.fq / g.cl
}


# Filter entries on cancer_max_rank
d.f <- d |>
  dplyr::filter(cancer_max_rank >= 0.5)


# Format table data
d.p <- d.f |>
  dplyr::rename(tsgene=tumor_suppressor) |>
  dplyr::mutate(hgnc_id=stringr::str_c('HGNC:', hgnc_id))


# Write to disk
readr::write_tsv(d.p, file=fs::path(PREFIX, 'prepared.tsv'))
