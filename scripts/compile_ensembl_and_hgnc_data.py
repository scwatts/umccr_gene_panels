#!/usr/bin/env python3
import argparse
import pathlib
import re
import sys


RE_ATTRIBUTE_GROUPS = re.compile(r'^(?P<name>[^"]+)(?: "(?P<data>.+)")?;?$')


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
    parser.add_argument('--annotations_fp', required=True, type=pathlib.Path)
    parser.add_argument('--output_dir', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.annotations_fp.exists():
        parser.error(f'Input file {args.annotations_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Gather Ensembl data
    ensembl_data = get_ensembl_data(args.annotations_fp)

    # Write as required
    write_gene_data(ensembl_data['genes'], ensembl_data['canonical_transcripts'], args.output_dir)
    write_cds_data(ensembl_data['cds'], ensembl_data['canonical_transcripts'], args.output_dir)


def get_ensembl_data(fp):
    genes = dict()
    cds = dict()
    canonical_transcripts = dict()
    with fp.open('r') as fh:
        # Skip header rows
        start_pos = int()
        while (line := fh.readline()):
            if not line.startswith('#'):
                break
            start_pos = fh.tell()
        fh.seek(start_pos)

        # Iterate GTF records; full pass prior to write required to map gene ID -> canonical transcript ID
        for i, line in enumerate(fh, 1):
            # Get record and skip entries without a HGNC ID
            record = prepare_record(line)
            if not record.hgnc_id:
                continue

            # Set in mapping with HGNC ID as key
            if record.feature == 'gene':

                # NOTE(SW): some genes are duplicated with the same HGNC ID and symbol. These duplicates are seen to
                # either have different Ensembl gene IDs with slightly different locations (e.g. TBCE [18 total]) or the
                # same Ensembl gene ID with a location on chrX and another in the PAR chrY region (e.g. SHOX [26 in
                # total]). Notably, these duplicates includes genes on the somatic panels; core panel: PPP2R3B; core
                # /and/ final panel: CRLF2, P2RY8.

                assert record.ensembl_gene_id not in genes
                genes[record.ensembl_gene_id] = record

            elif record.feature == 'CDS' and 'Ensembl_canonical' in record.attribute_dict.get('tag', {}):
                if record.ensembl_gene_id not in cds:
                    cds[record.ensembl_gene_id] = list()
                cds[record.ensembl_gene_id].append(record)

                if record.ensembl_gene_id in canonical_transcripts:
                    assert canonical_transcripts[record.ensembl_gene_id] == record.attribute_dict['transcript_id']
                canonical_transcripts[record.ensembl_gene_id] = record.attribute_dict['transcript_id']


    return {'genes': genes, 'cds': cds, 'canonical_transcripts': canonical_transcripts}


def prepare_record(line):
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

    # Set convenience attributes
    record.symbol = record.attribute_dict.get('gene_name')
    record.hgnc_id = record.attribute_dict.get('hgnc_id')
    record.ensembl_gene_id = record.attribute_dict.get('gene_id')
    record.ensembl_transcript_id = record.attribute_dict.get('transcript_id')

    return record


def write_gene_data(genes, canonical_transcripts, output_dir):
    header_tokens = (
        'hgnc_id',
        'ensembl_gene_id',
        'symbol',
        'ensembl_transcript_id',
        'contig',
        'start',
        'end',
        'strand',
    )

    genes_fp = output_dir / 'ensembl_gene_data.tsv'
    with genes_fp.open('w') as fh:
        print(*header_tokens, sep='\t', file=fh)
        for record in genes.values():
            data = [
                record.hgnc_id,
                record.ensembl_gene_id,
                record.symbol,
                canonical_transcripts.get(record.ensembl_gene_id, ''),
                record.seqname,
                record.start,
                record.end,
                record.strand,
            ]
            print(*data, sep='\t', file=fh)


def write_cds_data(cds, canonical_transcripts, output_dir):
    cds_fp = output_dir / 'ensembl_cds_data.tsv'
    with cds_fp.open('w') as fh:
        for records in cds.values():
            assert len({r.ensembl_transcript_id for r in records}) == 1
            for record in records:
                name = ';'.join((
                    record.symbol,
                    record.hgnc_id,
                    record.ensembl_gene_id,
                    record.ensembl_transcript_id,
                ))

                data = [
                    record.seqname,
                    record.start,
                    record.end,
                    name,
                    '.',  # BED score column
                    record.strand,
                ]
                print(*data, sep='\t', file=fh)


if __name__ == '__main__':
    main()
