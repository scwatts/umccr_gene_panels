# OncoVar

```bash
wget -P data/ https://oncovar.tania.wang:5443/resource/download/Onco_genes_OncoVar_TCGA/TCGA.PanCancer.onco.genes.OncoVar.tsv.gz
wget -P data/ https://oncovar.tania.wang:5443/resource/download/Onco_genes_OncoVar_ICGC/ICGC.PanCancer.onco.genes.OncoVar.tsv.gz

parallel gzip -d ::: data/*gz
```
