#!/usr/bin/env python3
import argparse
import pathlib


CHROM_ORDER = [f'chr{e}' for e in [*range(1, 23), 'X', 'Y']]


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--ensembl_gene_data_fp', required=True, type=pathlib.Path)
    parser.add_argument('--refseq_gene_data_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.ensembl_gene_data_fp.exists():
        parser.error(f'Input file {args.ensembl_gene_data_fp} does not exist')
    if not args.refseq_gene_data_fp.exists():
        parser.error(f'Input file {args.refseq_gene_data_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in gene data and panel data
    ensembl_genes = get_gene_data(args.ensembl_gene_data_fp, 'ensembl_gene_id')
    refseq_genes = get_gene_data(args.refseq_gene_data_fp, 'ncbi_gene_id')
    panel_records = get_panel_records(args.panel_fp)

    # Get gene data to write out
    data = list()
    seen_records = {'gene_records': set(), 'gene_ids': set()}
    for panel_record in panel_records:

        # Get Ensembl gene record where available, otherwise get RefSeq/NCBI record
        if panel_record['ensembl_gene_id'] != 'NA':
            gene_id = panel_record['ensembl_gene_id']
            gene_record = ensembl_genes[gene_id]
        else:
            gene_id = panel_record['ncbi_gene_id']
            gene_record = refseq_genes[gene_id]

        # Sanity check
        assert gene_record.values() not in seen_records['gene_records']
        assert gene_id not in seen_records['gene_ids']
        seen_records['gene_records'].add(gene_record.values())
        seen_records['gene_ids'].add(gene_id)

        # Print data
        data.append((
            gene_record['contig'],
            # NOTE(SW): BED start position is 0-based
            int(gene_record['start']) - 1,
            gene_record['end'],
            gene_record['symbol'],
        ))


    # Sort data by chromosome and location
    for d in sorted(data, key=lambda k: (CHROM_ORDER.index(k[0]), k[1])):
        print(*d, sep='\t')


def get_panel_records(fp):
    panel_records = list()
    with fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')

        for line in fh:
            data = line.rstrip().split('\t')
            panel_record = {k: v for k, v in zip(header_tokens, data)}
            panel_records.append(panel_record)

    return panel_records


def get_gene_data(fp, key_name):
    gene_records = dict()
    with fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')
        for line in fh:
            data = line.rstrip().split('\t')
            gene_record = {k: v for k, v in zip(header_tokens, data)}
            gene_records[gene_record[key_name]] = gene_record
    return gene_records


if __name__ == '__main__':
    main()
