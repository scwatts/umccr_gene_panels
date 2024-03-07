# Hartwig Medical Foundation (HMF)

```bash
wget -P data/ https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz

tar \
  -xzvf data/hmf_dna_pipeline_resources.38_v5.33.tar.gz \
  -s '#^.*/#data/#' \
  ./common/DriverGenePanel.38.tsv \
  ./sv/known_fusion_data.38.csv
```
