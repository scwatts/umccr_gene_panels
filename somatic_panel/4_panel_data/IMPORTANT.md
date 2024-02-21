# Important

## `umccr_cancer_genes.hg38.ensembl107.sort.bed`

* only used to select variants in hypermutated samples
* ***no slop***, check if need to compensate boundaries in bolt

## `umccr_cancer_genes.hg38.coding.bed`

* only used in allele frequency calculation (prior to cancer report)
* using all APPRIS principal and alternative transcripts
  * somewhat conflicting with Hartwig wrt use of Ensembl canonical transcripts

## hmftools panel data

* genes with ambiguous (oncogene and tsgene) or absent roles handled by:
  * deferring to Hartwig configuration where available
  * fallback to tsgene in all other cases
* impact of this is minimal since gene role is only used in driver likelihood calculation
  * not used by curators
* best solution would be to have hmftools be able to support dual role
* otherwise consult/resolve with curators if driver likelihood calculation is to be used
