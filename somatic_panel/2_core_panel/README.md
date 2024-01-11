# Core panel

Construct the core panel and generate the change set

```bash
./scripts/run.R

diff -u core_panel.tsv core_panel.new.tsv > core_panel.changes.tsv
```

Provide `core_panel.changes.tsv` to the curation team for review then make relevant changes to the curation gene lists
with their feedback. Once changes have been completed, rename `core_panel.new.tsv` to `core_panel.tsv` and delete
`core_panel.changes.tsv`.
