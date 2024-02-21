# Important

## `predispose_genes.hg38.transcript.bed`

* used to subset germline small variant calls
* critically, variants identified here may trigger referral
* hence for relevant genes we need to include all locations where there may be a relevant variant
* using all APPRIS principal and alternative transcripts, though this still does not guarantee all curatable variants
  * not currently generated; deferring to `umccr_predisposition_genes.genes.bed` to avoid stated potential issue
