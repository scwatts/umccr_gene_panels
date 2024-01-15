# UMCCR Cancer Gene Panels

In progress workspace for updating UMCCR gene panels.

## Panels

* [Somatic panel](somatic_panel/)
* [Fusion panel](fusion_panel/)
* [Germline panel](germline_panel/)

## Generating panels

## Common resources

### Gene annotations

Ensembl 105 and the HGNC are used to reference genes in a consistent way, enabling comparison between different sources.
This data is obtained from the GENCODE v39 release annotation GTF and processed so that construction of panels are
readily reproducible offline.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
gzip -d resources/gencode.v39.annotation.gtf.gz

./scripts/compile_ensembl_and_hgnc_data.py \
  --annotations_fp resources/gencode.v39.annotation.gtf \
  --output_dir resources/
```

I've anecdotally observed that gene symbols from sources match with a higher rate to the latest HGNC release rather than
the HGNC release tied to Ensembl 105. The general process taken is to then match gene symbols against the latest HGNC
release and subsequently add the Ensembl 105 gene ID. Similar to above, the file containing the latest HGNC data is
stored for offline use.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-11-01.tsv
```

### Hartwig Ensembl data cache and Gene Utilities

* TODO: write brief description
  * some panel data is built from the Hartwig Ensembl data cache
    * converting UMCCR panels to hmftools-format, generate fusion data for hmftools
  * this requires Gene Utilities from hmftools


```bash
mkdir -p resources/hmftools_ensembl_data_cache/

java \
  -cp /Users/stephen/projects/umccr_specific_hmftools_resources/1_small_variants/2_ensembl_105/software/gene-utils_a156ed6.jar \
  com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache \
    -ensembl_user anonymous \
    -ensembl_db mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_105_38 \
    -ref_genome_version 38 \
    -log_debug \
    -output_dir resources/hmftools_ensembl_data_cache/
```
