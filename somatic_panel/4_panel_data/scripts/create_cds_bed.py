#!/usr/bin/env python3
import argparse
import pathlib


CHROM_ORDER = [f'chr{e}' for e in [*range(1, 23), 'X', 'Y']]


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--ensembl_cds_data_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.ensembl_cds_data_fp.exists():
        parser.error(f'Input file {args.ensembl_cds_data_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in gene data and panel data
    ensembl_cds = get_cds_data(args.ensembl_cds_data_fp)
    panel_records = get_panel_records(args.panel_fp)

    # Some genes do not have CDS records or are completely absent from Ensembl 105
    skip_genes = {
        'ENSG00000130600.20', # H19
        'ENSG00000259104.2',  # PTCSC3
        'ENSG00000270123.4',  # VTRNA2-1
        'ENSG00000270141.3',  # TERC
        'ENSG00000280757.1',  # DUX4L1
        'ENSG00000284182.1',  # MIR143
        'NA',                 # TRA, TRB, IGH, TGK, TGL (no Ensembl record)
    }

    # Get CDS data to write out
    data = list()
    seen_transcript_ids = set()
    for panel_record in panel_records:

        if panel_record['ensembl_gene_id'] in skip_genes:
            continue

        assert panel_record['ensembl_transcript_id'] in ensembl_cds
        assert panel_record['ensembl_transcript_id'] not in seen_transcript_ids

        seen_transcript_ids.add(panel_record['ensembl_transcript_id'])
        cds_records = ensembl_cds[panel_record['ensembl_transcript_id']]

        assert len({r['attr_dict']['ensembl_transcript_id'] for r in cds_records}) == 1
        assert len({r['attr_dict']['ensembl_gene_id'] for r in cds_records}) == 1

        data.extend(cds_records)

    # Sort then output CDS data
    for record in sorted(data, key=cds_record_sort_key):
        d = (
            record['seqname'],
            # NOTE(SW): BED start position is 0-based
            int(record['start']) - 1,
            record['end'],
            record['attr'],
            '.', # BED score field
            record['strand'],
        )
        print(*d, sep='\t')


def cds_record_sort_key(r):
    return (CHROM_ORDER.index(r['seqname']), int(r['start']))


def get_panel_records(fp):
    panel_records = list()
    with fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')

        for line in fh:
            data = line.rstrip().split('\t')
            panel_record = {k: v for k, v in zip(header_tokens, data)}
            panel_records.append(panel_record)

    return panel_records


def get_cds_data(fp):
    header = ('seqname', 'start', 'end', 'attr', 'strand')
    attr_fields = ('symbol', 'hgnc_id', 'ensembl_gene_id', 'ensembl_transcript_id')

    records = dict()
    with fp.open('r') as fh:
        for line in fh:
            data = line.rstrip().split('\t')

            record = {k: v for k, v in zip(header, data)}
            record['attr_dict'] = {k: v for k, v in zip(attr_fields, record['attr'].split(';'))}
            key = record['attr_dict']['ensembl_transcript_id']

            if key not in records:
                records[key] = list()
            records[key].append(record)

    return records


if __name__ == '__main__':
    main()
