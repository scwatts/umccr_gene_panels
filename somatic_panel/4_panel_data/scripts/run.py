#!/usr/bin/env python3
import argparse
import collections
import pathlib
import re
import sys


import pickle


RE_ATTRIBUTE_GROUPS = re.compile(r'^(?P<name>[^"]+)(?: "(?P<data>.+)")?;?$')
# NOTE(SW): version and include non-numeric characters e.g. ENSG00000124333.16_PAR_Y
RE_ENSEMBL_GENE_ID_VERSION = re.compile(r'\..+$')


class GtfRecord:

    base_fields = (
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute',
    )

    def __init__(self, data_fields):
        for field_name, field_data in zip(self.base_fields, data_fields):
            setattr(self, field_name, field_data)

        self.attribute_dict = dict()
        self.ensembl_gene_id = str()
        self.hgnc_id = str()


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--annotations_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.annotations_fp.exists():
        parser.error(f'Input file {args.annotations_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get annotation data
    genes_ann, cds_ann = get_annotation_data(args.annotations_fp)

    # Read in gene panel
    panel_records = get_panel_data(args.panel_fp)

    # Get CDS
    missing_known_genes = {
        'HGNC:5477',  # IGH
        'HGNC:5715',  # IGK
        'HGNC:5853',  # IGL
        'HGNC:12027', # TRA
        'HGNC:12155', # TRB
    }

    missing_known_cds = {
        'HGNC:3082',  # DUX4L1
        'HGNC:4713',  # H19
        'HGNC:11727', # TERC
        'HGNC:31530', # MIR143
        'HGNC:37054', # VTRNA2-1
        'HGNC:43959', # PTCSC3
    }

    for record in panel_records:

        if False:
            # Skip missing
            if record.hgnc_id not in genes_ann and record.hgnc_id in missing_known_genes:
                continue

            assert record.hgnc_id in genes_ann
            d = genes_ann[record.hgnc_id]

            assert record.ensembl_gene_id == d.ensembl_gene_id

            print(record.symbol, d.ensembl_gene_id, record.hgnc_id, d.seqname, d.start, d.end, d.strand, sep='\t')



        # Skip missing
        if record.hgnc_id not in cds_ann and record.hgnc_id in missing_known_cds | missing_known_genes:
            continue

        if record.hgnc_id not in cds_ann:
            print(record.hgnc_id, record.symbol)
            continue



        assert record.hgnc_id in cds_ann
        ds = cds_ann[record.hgnc_id]

        ensembl_gene_ids = {d.ensembl_gene_id for d in ds}
        assert len(ensembl_gene_ids) == 1
        [ensembl_gene_id] = ensembl_gene_ids
        assert record.ensembl_gene_id == ensembl_gene_id

        for d in ds:
            print(record.symbol, d.ensembl_gene_id, record.hgnc_id, d.seqname, d.start, d.end, d.strand, sep='\t')


def get_panel_data(fh):
    with fh.open('r') as fh:
        header_fields = fh.readline().rstrip().split('\t')
        GeneRecord = collections.namedtuple('GeneRecord', header_fields)
        records = [GeneRecord(*line.rstrip().split('\t')) for line in fh]
    return records


def get_annotation_data(fp):

    with open('genes.pickle', 'rb') as fh:
        genes = pickle.load(fh)
    #cds = None
    with open('cds.pickle', 'rb') as fh:
        cds = pickle.load(fh)

    return genes, cds


    # Read gene bounds and Ensembl canonical transcript CDS
    with fp.open('r') as fh:
        # Skip header rows
        start_pos = int()
        while (line := fh.readline()):
            if not line.startswith('#'):
                break
            start_pos = fh.tell()
        fh.seek(start_pos)

        # Process
        genes = dict()
        cds = dict()
        for i, line in enumerate(fh, 1):

            # Report progress
            if (i % 10000 == 0):
                print(i, file=sys.stderr)

            # Create record and mapping for attribute data; some wasted cycles for records that we discard
            record = GtfRecord(line.rstrip('\n').split('\t'))

            attr_dict = dict()
            for attr_str in record.attribute.split('; '):
                # Get attribute entries, for entries with no data just set name -> True in dict
                attr_re = RE_ATTRIBUTE_GROUPS.match(attr_str)
                attr_name = attr_re['name']
                attr_data = attr_re['data'] if attr_re['data'] else True

                # The ont and tag attribute entries can occur multiple times, handle here
                if attr_name in {'ont', 'tag'}:
                    if attr_name not in attr_dict:
                        attr_dict[attr_name] = set()
                    attr_dict[attr_name].add(attr_data)
                else:
                    assert attr_name not in attr_dict
                    attr_dict[attr_name] = attr_data

            record.attribute_dict = attr_dict

            # Set Ensembl gene ID without version and HGNC ID
            record.ensembl_gene_id = RE_ENSEMBL_GENE_ID_VERSION.sub('', record.attribute_dict.get('gene_id'))
            record.hgnc_id = record.attribute_dict.get('hgnc_id')

            # Set in mapping with HGNC ID as key
            if record.feature == 'gene':
                genes[record.hgnc_id] = record
            elif record.feature == 'CDS' and 'Ensembl_canonical' in record.attribute_dict.get('tag', {}):
                if record.hgnc_id not in cds:
                    cds[record.hgnc_id] = list()
                cds[record.hgnc_id].append(record)

    with open('genes.pickle', 'wb') as fh:
        pickle.dump(genes, fh)
    with open('cds.pickle', 'wb') as fh:
        pickle.dump(cds, fh)

    return genes, cds


if __name__ == '__main__':
    main()
