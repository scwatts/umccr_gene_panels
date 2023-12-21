# UMCCR Cancer Gene Panels

In progress workspace for updating UMCCR gene panels.

## Common resources

Ensembl 105 and the associated HGNC release are used to identify genes, enabling comparison between different sources
following reconciliation. This data is obtained from the BioMart MySQL database dump rather than using the API so that
it is readily reproducible offline.

```bash
temp_dir=$(mktemp -d ./tmp.XXXXXXXX)

mkdir -p resources/ ${temp_dir}/

wget -P ${temp_dir}/ https://ftp.ensembl.org/pub/release-105/mysql/ensembl_mart_105/hsapiens_gene_ensembl__gene__main.txt.gz

# hsapiens_gene_ensembl__gene__main
#wget -P ${temp_dir}/ https://ftp.ensembl.org/pub/release-105/mysql/ensembl_mart_105/ensembl_mart_105.sql.gz

{
  echo -e "symbol\tensembl_gene_id\tname\thgnc_id";
  gzip -cd ${temp_dir}/hsapiens_gene_ensembl__gene__main.txt.gz | \
    grep -v 'CHR_HSCHR\|CHR_HG' | \
    awk -F$'\t' '
      BEGIN { OFS="\t"}
      $3 == "protein_coding" && $6 == "HGNC Symbol" { print $8, $32, $11}
    ' | \
    sed 's/\([^\t]\+\) \[Source.*HGNC:\([0-9]\+\).*$/\1\tHGNC:\2/g' | \
    sort;
} > resources/ensembl_105_genes_with_hgnc.tsv

rm -r ${temp_dir}
unset temp_dir
```

As a fallback gene symbols provided by sources are checked against the most recent HGNC release.

```bash
wget -P resources/ https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-11-01.tsv
```
