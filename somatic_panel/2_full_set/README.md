# Select entries from source gene lists

```bash
./scripts/run.R

diff -u gene_list.tsv gene_list.new.tsv > gene_list.changes.tsv
```

Provide `gene_list.changes.tsv` to curation team to review and make changes to the curation gene lists. Once curation
team have made their changes, rename `gene_list.new.tsv` to `gene_list.tsv` then delete `gene_list.changes.tsv`.
