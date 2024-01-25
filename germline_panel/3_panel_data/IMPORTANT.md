# Important

## `predispose_genes.hg38.transcript.bed`

* used to subset germline small variant calls
* critically, variants identified here may trigger referral
* hence, for relevant genes we need to include all locations where there may be a relevant variant
* using the entire length of the gene as defined by Ensembl 105 where available otherwise RefSeq
  * previously all APPRIS principal and alternative transcripts
