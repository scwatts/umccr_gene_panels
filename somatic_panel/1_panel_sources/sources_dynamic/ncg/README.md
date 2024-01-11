# Network of Cancer Genes (NCG)

```bash
mkdir -p data/

curl > data/NCG_cancerdrivers_annotation_supporting_evidence.tsv \
  -X POST \
  -d 'downloadcancergenes=Download' \
  http://network-cancer-genes.org/download.php
```
