#!/usr/bin/env python3
import argparse
import gzip
import pathlib
import re
import sys


CHROM_ORDER = [f'chr{e}' for e in [*range(1, 23), 'X', 'Y', 'M']]


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

        # Generic fields for writing
        self.gene_name = str()
        self.hgnc_id = str()
        self.gene_id = str()
        self.transcript_id = str()


def get_arguments():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='subcommand', dest='subcommand')

    parser_ensembl = subparsers.add_parser('ensembl')
    parser_ensembl.add_argument('--annotations_fp', required=True, type=pathlib.Path)
    parser_ensembl.add_argument('--appris_fp', required=True, type=pathlib.Path)
    parser_ensembl.add_argument('--output_dir', required=True, type=pathlib.Path)

    parser_refseq = subparsers.add_parser('refseq')
    parser_refseq.add_argument('--annotations_fp', required=True, type=pathlib.Path)
    parser_refseq.add_argument('--appris_fp', required=True, type=pathlib.Path)
    parser_refseq.add_argument('--contig_mapping_fp', required=True, type=pathlib.Path)
    parser_refseq.add_argument('--output_dir', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.annotations_fp.exists():
        parser.error(f'Input file {args.annotations_fp} does not exist')
    if not args.appris_fp.exists():
        parser.error(f'Input file {args.appris_fp} does not exist')
    if not args.subcommand != 'refseq' and not args.contig_mapping_fp.exists():
        parser.error(f'Input file {args.contig_mapping_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get requested annotations
    if args.subcommand == 'ensembl':
        annotations = compile_ensembl_data(args)
    elif args.subcommand == 'refseq':
        annotations = compile_refseq_data(args)
    else:
        assert False

    # Write data
    write_gene_data(annotations['genes'], args.subcommand, args.output_dir)
    write_transcript_data(annotations['transcripts'], args.subcommand, args.output_dir)
    write_cds_data(annotations['cds'], args.subcommand, args.output_dir)


def compile_ensembl_data(args):
    # Read in APPRIS annotations and relevant annotations
    appris_data = read_appris_data(args.appris_fp)
    annotations = retrieve_relevant_annotations(appris_data, args.annotations_fp)

    # Configure record attributes for generic write functions; modified inplace
    for record in (r for rs in annotations.values() for r in rs):
        [record.hgnc_id] = record.attribute_dict.get('hgnc_id', ['NA'])
        [record.gene_id] = record.attribute_dict['gene_id']
        [record.gene_name] = record.attribute_dict['gene_name']
        [record.transcript_id] = record.attribute_dict.get('transcript_id', ['NA'])

    return annotations


def compile_refseq_data(args):
    # Read in APPRIS annotations, contig data, relevant annotations
    appris_data = read_appris_data(args.appris_fp)
    contig_data = get_refseq_contig_data(args.contig_mapping_fp)
    annotations = retrieve_relevant_annotations(appris_data, args.annotations_fp)

    # Similar to Ensembl annotations I want HGNC in associated with transcripts and CDS + stop
    # codons, however RefSeq only provides this data in gene features. Here I lazily iterate the
    # entire record set another time to collect all RefSeq HGNC tags so I can assign below
    hgnc_ids = dict()
    for record in (r for rs in annotations.values() for r in rs):
        if (hgnc_id := get_refseq_dbxref_data(record, 'HGNC')) == 'NA':
            continue

        # NOTE(SW): not assigning here for consistency, retrived again below
        gene_id = get_refseq_dbxref_data(record, 'GeneID')

        assert gene_id != 'NA'
        if gene_id in hgnc_ids:
            assert hgnc_ids[gene_id] == hgnc_id

        hgnc_ids[gene_id] = hgnc_id

    # Configure record attributes for generic write functions; modified inplace
    for record in (r for rs in annotations.values() for r in rs):
        if not (seqname := contig_data.get(record.seqname)):
            continue
        record.seqname = seqname

        record.gene_id = get_refseq_dbxref_data(record, 'GeneID')
        [record.gene_name] = record.attribute_dict['gene']
        record.transcript_id = record.attribute_dict.get('transcript_id')[0] or 'NA'

        if not record.hgnc_id:
            record.hgnc_id = hgnc_ids.get(record.gene_id, 'NA')

    # Another pass to exclude records that have non-main contigs
    contigs = set(contig_data.values())
    for key in annotations.keys():
        annotations[key] = [r for r in annotations[key] if r.seqname in contigs]

    return annotations


def read_appris_data(fp):
    records = dict()
    with fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')
        for line in fh:
            record = {k: v for k, v in zip(header_tokens, line.rstrip().split('\t'))}

            if not (annotation := record.get('APPRIS Annotation')):
                continue

            if annotation.startswith('ALTERNATIVE') or annotation.startswith('PRINCIPAL'):
                assert record['Transcript ID'] not in records
                records[record['Transcript ID']] = record
    return records


def get_refseq_contig_data(fp):
    data = dict()
    with fp.open('r') as fh:
        for line in fh:
            name, accession = line.rstrip().split('\t')
            assert accession not in data
            data[accession] = name
    return data


def retrieve_relevant_annotations(appris_data, fp):
    genes = list()
    transcripts = list()
    cds = list()
    for i, record in enumerate(gtf_record_iterator(fp), 1):

        if i % 10000 == 0:
            print(i, file=sys.stderr)

        # Get relevant features
        # NOTE(SW): start codons always fall within CDS; checked manually
        if record.feature == 'gene':
            genes.append(record)
        elif record.feature in {'transcript', 'CDS', 'stop_codon'}:

            # Select only APPRIS transcripts
            if not (transcript_full_id := record.attribute_dict.get('transcript_id')):
                continue

            assert len(transcript_full_id) == 1
            [transcript_full_id] = transcript_full_id
            transcript_id = re.sub(r'\.\d+$', '', transcript_full_id)

            # NOTE(SW): APPRIS has versioned transcripts for RefSeq but not Ensembl, so I check both
            if transcript_id not in appris_data and transcript_full_id not in appris_data:
                continue

            if record.feature == 'transcript':
                transcripts.append(record)
            elif record.feature in {'CDS', 'stop_codon'}:
                cds.append(record)
            else:
                assert False

    return {'genes': genes, 'transcripts': transcripts, 'cds': cds}


def gtf_record_iterator(fp):
    with gzip.open(fp, 'rt') as fh:
        # Skip header rows
        start_pos = int()
        while (line := fh.readline()):
            if not line.startswith('#'):
                break
            start_pos = fh.tell()
        fh.seek(start_pos)

        for line in fh:

            # NOTE(SW): the final line in the RefSeq GTF is '###', must handle here...
            if line == '###\n':
                continue

            yield prepare_record(line)


def prepare_record(line):
    # Create record and mapping for attribute data
    # NOTE(SW): some wasted cycles for records that we discard
    record = GtfRecord(line.rstrip('\n').split('\t'))

    attr_dict = dict()
    for attr_name, attr_data in attr_string_parser(record.attribute):

        # NOTE(SW): use list types for all attribute values to implicitly handle duplicate names
        if attr_name not in attr_dict:
            attr_dict[attr_name] = list()
        attr_dict[attr_name].append(attr_data)

    record.attribute_dict = attr_dict

    return record


def attr_string_parser(s):
    # NOTE(SW): custom token parser required to cover all following cases:
    #   * quoted data
    #   * unquoted data
    #   * special characters (double quote, semicolon, spaces) in quoted data
    #   * spaces in unquoted data

    # This, unsurprisingly, is painfully slow in Python

    attr_name = list()
    attr_data = list()

    name_flag = True
    quote_flag = False
    escape_flag = False
    end_flag = False

    for l in s:

        if l == ';' and not quote_flag:
            yield ''.join(attr_name), ''.join(attr_data)
            attr_name = list()
            attr_data = list()
            name_flag = True
            end_flag = True
        elif l == ' ':
            if end_flag:
                end_flag = False
            elif name_flag:
                name_flag = False
            else:
                attr_data.append(l)
        elif l == '"':
            if escape_flag:
                attr_data.append(l)
                escape_flag = False
            else:
                quote_flag = not quote_flag
        elif l == '\\':
            escape_flag = True
        elif name_flag:
            attr_name.append(l)
        else:
            attr_data.append(l)

        if escape_flag and l != '\\':
            escape_flag = False

    return attr_name, attr_data


def get_refseq_dbxref_data(record, name):
    data = [t for t in record.attribute_dict['db_xref'] if t.startswith(name)]
    assert 0 <= len(data) <= 1
    return data[0].replace(f'{name}:', '', 1) if data else 'NA'


def write_gene_data(data, source, output_dir):
    header = ['hgnc_id', 'gene_id', 'symbol', 'contig', 'start', 'end', 'strand']
    with (output_dir / f'{source}.genes.tsv').open('w') as fh:
        print(*header, sep='\t', file=fh)
        for record in sorted(data, key=gtf_record_sort):
            print(
                record.hgnc_id,
                record.gene_id,
                record.gene_name,
                record.seqname,
                record.start,
                record.end,
                record.strand,
                sep='\t',
                file=fh,
            )

def write_transcript_data(data, source, output_dir):
    with (output_dir / f'{source}.transcripts.bed').open('w') as fh:
        for record in sorted(data, key=gtf_record_sort):

            hgnc_id = '' if record.hgnc_id == 'NA' else record.hgnc_id
            info_fields = [record.gene_name, hgnc_id, record.gene_id, record.transcript_id]

            print(
                record.seqname,
                int(record.start)-1,
                record.end,
                ';'.join(info_fields),
                record.strand,
                sep='\t',
                file=fh,
            )


def write_cds_data(data, source, output_dir):
    with (output_dir / f'{source}.cds.bed').open('w') as fh:
        for record in sorted(data, key=gtf_record_sort):

            hgnc_id = '' if record.hgnc_id == 'NA' else record.hgnc_id
            info_fields = [record.gene_name, hgnc_id, record.gene_id, record.transcript_id, record.feature]

            print(
                record.seqname,
                int(record.start)-1,
                record.end,
                ';'.join(info_fields),
                record.strand,
                sep='\t',
                file=fh,
            )


def gtf_record_sort(r):
    return (CHROM_ORDER.index(r.seqname), int(r.start))


if __name__ == '__main__':
    main()
