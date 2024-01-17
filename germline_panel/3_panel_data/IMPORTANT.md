# Important

## `predispose_genes.hg38.transcript.bed`

* used to subset germline small variant calls
* critically, variants identified here may trigger referral
* hence, for relevant genes we need to include all locations where there may be a relevant variant
* imo this should be the entire length of the gene, not just CDS
* currently using only Ensembl 105 canonical transcripts
* not including stop codon in regions
