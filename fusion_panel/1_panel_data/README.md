# Panel data

Get Hartwig fusion data

```bash
wget -P data/ https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz

tar \
  -xzvf data/hmf_dna_pipeline_resources.38_v5.33.tar.gz \
  -s '#^.*/#data/#' \
  ./sv/known_fusion_data.38.csv
```

Modify to align with UMCCR fusion panel and create input file for gene-utils. Two genes are not includes:

* AL121790.1 has no entry in HGNC but does correspond to ENSG00000258414 in the GRCh37 Ensembl 105 release. Hartwig have
  taken coordinates from that Ensembl entry and lifted over to GRCh38 in order to manually include this gene. See
  [here](https://github.com/hartwigmedical/hmftools/blob/master/gene-utils/README.md#global-gene-panel)
* C19MC also has no entry in HGNC, nor in Ensembl


```bash
mkdir -p output/

./scripts/run.R
```

Regenerate the Hartwig fusion data file with UMCCR gene entries and Ensembl 105 for use with sash or oncoanalyser

```bash
# NOTE(SW): hmftools Ensembl data cache expected to be in specific directory structure
db=$(mktemp -d tmp.XXXXXXXXXX)
dn=${db}/ensembl_data_cache/38
mkdir -p ${dn}/

ln -s $(find $(pwd)/../../resources/hmftools_ensembl_data_cache/ -type f) ${dn}/

# NOTE(SW): these genes are not available in the Ensembl 105 data cache
sed \
  -i \
  -e '/PROMISCUOUS_3.*\tIGH\t/d' \
  -e '/PROMISCUOUS_3.*\tIGK\t/d' \
  -e '/PROMISCUOUS_3.*\tIGL\t/d' \
  -e '/PROMISCUOUS_[35].*\tTRA\t/d' \
  -e '/PROMISCUOUS_[35].*\tTRB\t/d' \
  output/fusion_database.tsv

java \
  -cp ../../other/gene-utils-1.1-jar-with-dependencies.jar \
  com.hartwig.hmftools.geneutils.fusion.GenerateFusionFiles \
    -known_fusion_db_file output/fusion_database.tsv \
    -resource_repo_dir ${db}/ \
    -log_debug \
    -output_dir output/

rm -r ${db}/
```

Create fusion panel data required by `prioritize_sv.py`

```bash
{
  echo five_gene,three_gene;
  awk -F, 'BEGIN { OFS="," }; NR > 1 && ($2 && $3) { print $2, $3 }' output/hmftools_fusion_data.subset.csv | \
    sort | \
    uniq;
} > output/hartwig_known_pairs.csv

{
  echo gene;
  grep PROMISCUOUS output/hmftools_fusion_data.subset.csv | \
    awk -F, '$1 ~ /_5$|^IG_/ { print $2 }' | \
    sort | \
    uniq;
} > output/hartwig_promiscuous_five.csv

{
  echo gene;
  grep PROMISCUOUS output/hmftools_fusion_data.subset.csv | \
    awk -F, '$1 ~ /_3$|__TARGET$/ { print $3 }' | \
    sort | \
    uniq;
} > output/hartwig_promiscuous_three.csv
```
