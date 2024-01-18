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

./scripts/compile_ensembl_and_hgnc_data.py \
  --annotations_fp resources/gencode.v39.annotation.gtf.gz \
  --output_dir resources/
```

I've anecdotally observed that gene symbols from sources match with a higher rate to the latest HGNC release rather than
the HGNC release tied to Ensembl 105. The general process taken is to then match gene symbols against the latest HGNC
release and subsequently add the Ensembl 105 gene ID. Similar to above, the file containing the latest HGNC data is
stored for offline use.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-11-01.tsv
```

Lastly, some genes are not present in Ensembl 105 (e.g. IGH, TRA) and to obtain genomic coordinates I use RefSeq
annotations. Similar approach as Ensembl above.

```bash
base_ftp=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14

wget -P resources/ ${base_ftp}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
curl ${base_ftp}/GCF_000001405.40_GRCh38.p14_assembly_report.txt | \
  awk 'BEGIN { OFS="\t" } $0 ~ /^MT|^[0-9XY]/ { print "chr" $1, $7 }' > resources/refseq_contig_id_mapping.tsv

curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.summary.txt.gz | \
  gzip -cd | \
  awk -F$'\t' 'BEGIN { OFS="\t" } NR > 1 && $10 == "MANE Select" { print $1, $6 }' | \
  sed 's/^GeneID://' > resources/refseq_mane_select.tsv

./scripts/compile_refseq_data.py \
  --annotations_fp resources/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
  --contig_mapping_fp resources/refseq_contig_id_mapping.tsv \
  --mane_select_fp resources/refseq_mane_select.tsv \
  --output_dir resources/
```

### Hartwig Ensembl data cache and Gene Utilities

We use the hmftools workflows or components from Hartwig in our post-processing pipelines, and this requires specific
panel data files to be generated from the UMCCR gene lists. Generating the these files is done using `gene-utils` and is
performed for the somatic panel and fusion panel. Beyond the panel data files, we must also build the hmftools Ensembl
data cache for both the hmftools workflows and certain `gene-utils` functions.

First compile `gene-utils` at commit `a156ed6` (current v1.0 release is missing important changes), see
[COMPILE_GENE_UTILS.md](COMPILE_GENE_UTILS.md).

Build the Ensembl data cache

```bash
mkdir -p resources/hmftools_ensembl_data_cache/

java \
  -cp other/gene-utils_a156ed6.jar \
  com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache \
    -ensembl_user anonymous \
    -ensembl_db mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_105_38 \
    -ref_genome_version 38 \
    -log_debug \
    -output_dir resources/hmftools_ensembl_data_cache/
```
