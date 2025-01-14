# Resource creation

## Gene annotations

First download APPRIS annotations so that we can select the most relevant transcripts.

```bash
base_url=https://apprisws.bioinfo.cnio.es/pub/releases/2023_08.v48/datafiles/homo_sapiens

wget -O resources/appris.e105v46.tsv ${base_url}/e105v46/appris_data.appris.txt
wget -O resources/appris.rs110v48.tsv ${base_url}/rs110v48/appris_data.appris.txt
```

Annotations from Ensembl 105, HGNC, and RefSeq are used to reference genes in a consistent way, enabling exact and
confident comparison between different sources. The Ensembl 105 annotations are obtained from the GENCODE v39 release
and processed to be readily usable.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

./scripts/compile_annotation_data.py \
  ensembl \
  --annotations_fp resources/gencode.v39.annotation.gtf.gz \
  --appris_fp resources/appris.e105v46.tsv \
  --output_dir ./resources/
```

As some genes are not present in Ensembl 105 (e.g. IGH, TRA), I use RefSeq v110 annotation release inplace of Ensembl
when necessary. As RefSeq annotations use accessions for contigs, these must be renamed with using a lookup table.

```bash
base_ftp=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/GCF_000001405.40_GRCh38.p14
wget -P resources/ ${base_ftp}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

curl ${base_ftp}/GCF_000001405.40_GRCh38.p14_assembly_report.txt | \
  awk 'BEGIN { OFS="\t" } $0 ~ /^MT|^[0-9XY]/ { print "chr" $1, $7 }' | \
  sed 's/chrMT/chrM/' > resources/refseq_contig_id_mapping.tsv

./scripts/compile_annotation_data.py \
  refseq \
  --annotations_fp resources/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
  --appris_fp resources/appris.rs110v48.tsv \
  --contig_mapping_fp resources/refseq_contig_id_mapping.tsv \
  --output_dir ./resources/
```

I've anecdotally observed that source gene symbols match at a higher rate to the latest HGNC database release compared
to the release tied with Ensembl 105. The general process taken is to match gene symbols against the latest HGNC
database release and subsequently add the Ensembl 105 gene annotations. Hence, we retrive the latest HGNC complete set.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-11-01.tsv
```

## Hartwig Ensembl data cache and Gene Utilities

We use hmftools workflows/components from Hartwig in our post-processing pipelines, and this requires specific data
files to be generated from the UMCCR gene panels. These data files are created with `gene-utils` for the somatic panel
and fusion panel. Beyond the panel data files, we must also build the hmftools Ensembl data cache for both the hmftools
workflows and certain `gene-utils` functions.

First compile `gene-utils` at commit `a156ed6` (current v1.0 release is missing important changes), see
[COMPILE_GENE_UTILS.md](COMPILE_GENE_UTILS.md).

Build the Ensembl data cache for GRCh38 from the Ensembl 105 release

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
