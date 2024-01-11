# Panel data

* coding BED
* gene location BED

`../../resources/gencode.v39.annotation.gtf`
`tag "Ensembl_canonical";`


```bash
./scripts/run.py \
  --panel_fp ../3_final_panel/final_panel.tsv \
  --annotations_fp ../../resources/gencode.v39.annotation.gtf
```
