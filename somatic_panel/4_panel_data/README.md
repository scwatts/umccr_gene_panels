# Panel data

Generate panel data for UMCCR post-processing

```bash
mkdir -p output/

awk -F$'\t' 'NR > 1 { print ($1 != "NA" ? $1 : $3)  }' ../3_final_panel/final_panel.tsv > output/umccr_cancer_genes.latest.genes
awk -F$'\t' '$8 == "TRUE" { print ($1 != "NA" ? $1 : $3) }' ../3_final_panel/final_panel.tsv > output/umccr_cancer_genes.tsgenes.latest.genes

./scripts/create_gene_bed.py > output/umccr_cancer_genes.genes.bed \
  --panel_fp ../3_final_panel/final_panel.tsv \
  --ensembl_gene_data_fp ../../resources/ensembl_gene_data.tsv \
  --refseq_gene_data_fp ../../resources/refseq_gene_data.tsv

./scripts/create_cds_bed.py > output/umccr_cancer_genes.cds.bed \
  --panel_fp ../3_final_panel/final_panel.tsv \
  --ensembl_cds_data_fp ../../resources/ensembl.cds.bed \
  --refseq_cds_data_fp ../../resources/refseq.cds.bed
```

Generate hmftools-compatible panel data

```bash
wget -P data/ https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz

tar \
  -xzvf data/hmf_dna_pipeline_resources.38_v5.33.tar.gz \
  -s '#^.*/#data/#' \
  ./variants/clinvar.38.vcf.gz

db=$(mktemp -d tmp.XXXXXXXXXX)
mkdir -p ${db}/

./scripts/create_hmftools_panel.py > ${db}/umccr_cancer_genes.driver_panel.tsv \
  --panel_fp ../3_final_panel/final_panel.tsv \
  --hartwig_panel_fp ../1_panel_sources/sources_dynamic/hmf/data/DriverGenePanel.38.tsv


mkdir -p ${db}/resources/{ensembl_data_cache,sage}/38/

ln -s $(pwd)/../../resources/hmftools_ensembl_data_cache/* ${db}/resources/ensembl_data_cache/38/
ln -s $(pwd)/data/clinvar.38.vcf.gz ${db}/resources/sage/38/

java \
  -cp ../../other/gene-utils-1.1-jar-with-dependencies.jar \
  com.hartwig.hmftools.geneutils.drivers.GenerateDriverGeneFiles \
    -resource_repo_dir ${db}/resources/ \
    -driver_gene_panel ${db}/umccr_cancer_genes.driver_panel.tsv \
    -log_debug \
    -output_dir ${db}/output/

mkdir -p output/hmftools/

mv $(find ${db}/output/ -type f) output/hmftools/

rm -r ${db}/
```
