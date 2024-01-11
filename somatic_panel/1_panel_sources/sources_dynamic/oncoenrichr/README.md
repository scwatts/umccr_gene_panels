# oncoEnrichR

```bash
wget -P data/ http://insilico.hpc.uio.no/oncoEnrichR/v1.4.2/genedb_v1.4.2.rds

R --vanilla <<EOF
library(readr)
d <- readRDS('data/genedb_v1.4.2.rds')
readr::write_tsv(d\$all, 'data/genedb_v1.4.2.tsv')
EOF
```
