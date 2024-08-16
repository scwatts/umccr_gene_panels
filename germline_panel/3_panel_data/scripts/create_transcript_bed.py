#!/usr/bin/env python3
import argparse
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--ensembl_transcript_data_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.ensembl_transcript_data_fp.exists():
        parser.error(f'Input file {args.ensembl_transcript_data_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read gene IDs for panel data
    gene_ids = set()
    with args.panel_fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')

        for line in fh:
            data = line.rstrip().split('\t')
            panel_record = {k: v for k, v in zip(header_tokens, data)}
            gene_ids.add(panel_record['ensembl_gene_id'])

    # Store transcript to print later
    transcripts = list()
    genes_seen = set()
    header_tokens = ('chrom', 'start', 'stop', 'data_str', 'strand')
    data_tokens = ('symbol', 'hgnc_id', 'gene_id', 'transcript_id')
    with args.ensembl_transcript_data_fp.open('r') as fh:
        for line in fh:
            values = line.rstrip().split('\t')
            record = {k: v for k, v in zip(header_tokens, values)}
            record['data'] = {k: v for k, v in zip(data_tokens, record['data_str'].split(';'))}

            if record['data']['gene_id'] not in gene_ids:
                continue

            genes_seen.add(record['data']['gene_id'])
            transcripts.append(line.rstrip())

    # Check we have processed all genes
    genes_missed = gene_ids - genes_seen
    if genes_missed:
        plurality = 'genes' if len(genes_missed) > 1 else 'gene'
        genes_missed_strs = [f'  - {n}' for n in genes_missed]
        message = f'ERROR: no transcripts found for the following {plurality}:'
        print(message, *genes_missed_strs, sep='\n', file=sys.stderr)
        sys.exit(1)

    # Print out transcripts
    print(*transcripts, sep='\n')


if __name__ == '__main__':
    main()
