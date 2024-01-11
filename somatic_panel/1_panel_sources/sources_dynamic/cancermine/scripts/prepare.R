library(dplyr)
library(fs)
library(ggplot2)
library(readr)


# Read in table
d <- readr::read_delim(fs::path(PREFIX, 'data', 'cancermine_collated.tsv'), col_types='cccccccd')


# Review citation count distribution to empirically decide threshold; using lenient threshold of >=5 citations
if (FALSE) {
  library(ggplot2)
  library(patchwork)

  threshold <- 5

  g.fq <- ggplot2::ggplot(d[d$citation_count<=200, ], ggplot2::aes(x=citation_count))
  g.fq <- g.fq + ggplot2::geom_freqpoly(binwidth=1)
  g.fq <- g.fq + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.fq <- g.fq + ggplot2::xlim(1, NA)
  g.fq <- g.fq + ggplot2::theme_bw()
  g.fq <- g.fq + ggplot2::theme(axis.title.x=ggplot2::element_blank())
  g.fq <- g.fq + ggplot2::ggtitle('Frequency distribution')


  g.cl <- ggplot2::ggplot(d[d$citation_count<=200, ], ggplot2::aes(x=citation_count))
  g.cl <- g.cl + ggplot2::stat_ecdf(geom='step')
  g.cl <- g.cl + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.cl <- g.cl + ggplot2::theme_bw()
  g.cl <- g.cl + ggplot2::ggtitle('Cumulative distribution')
  g.cl <- g.cl + ggplot2::labs(x='Citations')

  g.fq / g.cl
}


# Exclude 'Driver' entries and filter records on citation count
d.f <- d |>
  dplyr::filter(role != 'Driver') |>
  dplyr::filter(citation_count >= 5)


# Process records for each gene
d.p <- d.f |>
  dplyr::group_by(gene_hugo_id) |>
  dplyr::summarise(
    oncogene='Oncogene' %in% role,
    tsgene='Tumor_Suppressor' %in% role,
    citation_count=sum(citation_count),
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(hgnc_id=gene_hugo_id)


# Review citation counts for second pass filter; using threshold of >=8 citations
if (FALSE) {
  threshold <- 8

  g.fq <- ggplot2::ggplot(d.p[d.p$citation_count<=50, ], ggplot2::aes(x=citation_count))
  g.fq <- g.fq + ggplot2::geom_freqpoly(binwidth=1)
  g.fq <- g.fq + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.fq <- g.fq + ggplot2::xlim(5, NA)
  g.fq <- g.fq + ggplot2::theme_bw()
  g.fq <- g.fq + ggplot2::theme(axis.title.x=ggplot2::element_blank())
  g.fq <- g.fq + ggplot2::ggtitle('Frequency distribution')


  g.cl <- ggplot2::ggplot(d.p[d.p$citation_count<=50, ], ggplot2::aes(x=citation_count))
  g.cl <- g.cl + ggplot2::stat_ecdf(geom='step')
  g.cl <- g.cl + ggplot2::geom_vline(xintercept=threshold, colour='red')
  g.cl <- g.cl + ggplot2::theme_bw()
  g.cl <- g.cl + ggplot2::ggtitle('Cumulative distribution')
  g.cl <- g.cl + ggplot2::labs(x='Citations')

  g.fq / g.cl
}


# Filter records on citation count
d.d <- d.p |> dplyr::filter(citation_count >= 8)


# Write to disk
readr::write_tsv(d.d, file=fs::path(PREFIX, 'prepared.tsv'))
