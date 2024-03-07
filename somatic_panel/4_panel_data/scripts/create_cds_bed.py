#!/usr/bin/env python3
import argparse
import pathlib


CHROM_ORDER = [f'chr{e}' for e in [*range(1, 23), 'X', 'Y']]


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--ensembl_cds_data_fp', required=True, type=pathlib.Path)
    parser.add_argument('--refseq_cds_data_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.ensembl_cds_data_fp.exists():
        parser.error(f'Input file {args.ensembl_cds_data_fp} does not exist')
    if not args.refseq_cds_data_fp.exists():
        parser.error(f'Input file {args.refseq_cds_data_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in gene data and panel data
    ensembl_cds = get_cds_data(args.ensembl_cds_data_fp)
    refseq_cds = get_cds_data(args.refseq_cds_data_fp)
    panel_records = get_panel_records(args.panel_fp)

    # Some genes do not have CDS records or are completely absent from Ensembl 105 and RefSeq 110
    skip_genes = {
        'HGNC:3082',  # DUX4L1
        'HGNC:4713',  # H19
        'HGNC:5477',  # IGH
        'HGNC:5715',  # IGK
        'HGNC:5853',  # IGL
        'HGNC:9709',  # PVT1
        'HGNC:11727', # TERC
        'HGNC:12027', # TRA
        'HGNC:12155', # TRB
        'HGNC:31530', # MIR143
        'HGNC:37054', # VTRNA2-1
        'HGNC:43959', # PTCSC3
    }

    # Get CDS data to write out
    data = list()
    seen_ids = set()
    for panel_record in panel_records:

        if panel_record['hgnc_id'] in skip_genes:
            assert panel_record['hgnc_id'] not in ensembl_cds
            assert panel_record['hgnc_id'] not in refseq_cds
            continue

        assert panel_record['hgnc_id'] not in seen_ids
        seen_ids.add(panel_record['hgnc_id'])

        if panel_record['hgnc_id'] in ensembl_cds:
            cds_records = ensembl_cds[panel_record['hgnc_id']]
        elif panel_record['hgnc_id'] in refseq_cds:
            # NOTE(SW): there are currently no entries where RefSeq but not Ensembl has CDS annotations
            cds_records = refseq_cds[panel_record['hgnc_id']]
        else:
            assert False

        assert len({r['attr_dict']['hgnc_id'] for r in cds_records}) == 1
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
    attr_fields = ('symbol', 'hgnc_id', 'gene_id', 'transcript_id', 'feature_type')

    records = dict()
    with fp.open('r') as fh:
        for line in fh:
            data = line.rstrip().split('\t')

            record = {k: v for k, v in zip(header, data)}
            record['attr_dict'] = {k: v for k, v in zip(attr_fields, record['attr'].split(';'))}
            key = record['attr_dict']['hgnc_id']

            if key not in records:
                records[key] = list()
            records[key].append(record)

    return records


if __name__ == '__main__':
    main()
