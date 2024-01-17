#!/usr/bin/env python3
import argparse
import gzip
import pathlib
import re


RE_HGNC_ID = re.compile(r'db_xref "HGNC:(HGNC:[0-9]+)";')
RE_GENE_SYMBOL = re.compile(r'gene_id "([^"]+)";')
RE_NCBI_GENE_ID = re.compile(r'db_xref "GeneID:([0-9]+)";')
RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)";')


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

            self.hgnc_id = str()
            self.symbol = str()
            self.ncbi_gene_id = str()
            self.mane_transcript_id = str()
            self.chrom = str()


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotations_fp', required=True, type=pathlib.Path)
    parser.add_argument('--contig_mapping_fp', required=True, type=pathlib.Path)
    parser.add_argument('--mane_select_fp', required=True, type=pathlib.Path)
    parser.add_argument('--output_dir', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.annotations_fp.exists():
        parser.error(f'Input file {args.annotations_fp} does not exist')
    if not args.contig_mapping_fp.exists():
        parser.error(f'Input file {args.contig_mapping_fp} does not exist')
    if not args.mane_select_fp.exists():
        parser.error(f'Input file {args.mane_select_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Get data
    contig_data = get_contig_data(args.contig_mapping_fp)
    mane_select_data = get_mane_select_data(args.mane_select_fp)
    refseq_data = get_refseq_data(args.annotations_fp, contig_data)

    # Write as required
    write_gene_data(refseq_data['genes'], mane_select_data, args.output_dir)
    # NOTE(SW): all transcripts are currently covered by Ensembl
    #write_cds_data(refseq_data['cds'], refseq_data['genes'], mane_select_data, args.output_dir)


def get_refseq_data(fp, contig_data):
    genes = dict()
    cds = dict()
    with gzip.open(fp, 'rt') as fh:
        for line in fh:
            # NOTE(SW): comments lines at start and end of file, just handle here in every
            # iteration of each line
            if line.startswith('#'):
                continue

            record = GtfRecord(line.rstrip('\n').split('\t'))

            # Lazily use regex to pull out relevant attribute data
            record.hgnc_id = get_attr_data(RE_HGNC_ID, record.attribute)
            record.ncbi_gene_id = get_attr_data(RE_NCBI_GENE_ID, record.attribute)
            record.symbol = get_attr_data(RE_GENE_SYMBOL, record.attribute)
            record.refseq_transcript_id = get_attr_data(RE_TRANSCRIPT_ID, record.attribute)

            # Convert chromosome accession to name, skip non-main records
            if record.seqname not in contig_data:
                continue
            record.chrom = contig_data[record.seqname]

            if record.feature == 'gene':

                assert record.ncbi_gene_id
                assert record.ncbi_gene_id != 'NA'

                # Exclude chrY PAR genes; transcripts don't need exclusion since use of MANE select
                # implicitly does that
                # NOTE(SW): manually ensured that qualifying via the '_1' suffix and chrY location
                # is a reliable and safe method to identify such genes
                if record.ncbi_gene_id in genes:
                    if record.chrom == 'chrY' and record.symbol.endswith('_1'):
                        continue
                    else:
                        assert False

                genes[record.ncbi_gene_id] = record

            elif record.feature == 'CDS':

                assert record.refseq_transcript_id
                assert record.refseq_transcript_id != 'NA'

                if record.refseq_transcript_id not in cds:
                    cds[record.refseq_transcript_id] = list()
                cds[record.refseq_transcript_id].append(record)

    return {'genes': genes, 'cds': cds}


def get_contig_data(fp):
    data = dict()
    with fp.open('r') as fh:
        for line in fh:
            name, accession = line.rstrip().split('\t')
            assert accession not in data
            data[accession] = name
    return data


def get_mane_select_data(fp):
    data = dict()
    with fp.open('r') as fh:
        for line in fh:
            ncbi_gene_id, refseq_transcript_id = line.rstrip().split('\t')
            assert ncbi_gene_id not in data
            data[ncbi_gene_id] = refseq_transcript_id
    return data


def get_attr_data(regex, attr_str):
    if (re_result := regex.search(attr_str)):
        return re_result.group(1)
    else:
        return 'NA'


def write_gene_data(genes, mane_transcripts, output_dir):
    header_tokens = (
        'hgnc_id',
        'ncbi_gene_id',
        'symbol',
        'mane_transcript_id',
        'contig',
        'start',
        'end',
        'strand',
    )

    genes_fp = output_dir / 'refseq_gene_data.tsv'
    with genes_fp.open('w') as fh:
        print(*header_tokens, sep='\t', file=fh)
        for record in genes.values():
            data = [
                record.hgnc_id,
                record.ncbi_gene_id,
                record.symbol,
                mane_transcripts.get(record.ncbi_gene_id, ''),
                record.chrom,
                record.start,
                record.end,
                record.strand,
            ]
            print(*data, sep='\t', file=fh)


def write_cds_data(cds, genes, mane_transcripts, output_dir):
    cds_fp = output_dir / 'mane_cds_data.tsv'
    with cds_fp.open('w') as fh:

        for gene_record in genes.values():
            # Some gene entries don't have transcripts
            if gene_record.ncbi_gene_id not in mane_transcripts:
                continue

            # Grab MANE select transcript records for relevant genes
            mane_transcript_id = mane_transcripts[gene_record.ncbi_gene_id]

            # Skipping MANE select transcripts that don't have CDS, this is caused by the
            # transcript being on a non-main contig
            # NOTE(SW): not fully explored for other corner cases
            if mane_transcript_id not in cds:
                continue
            transcript_records = cds[mane_transcript_id]

            # Ensure all CDS records are for the correct gene/transcript
            ncbi_gene_id = str()
            mane_transcript_id = str()
            for record in transcript_records:
                if not ncbi_gene_id:
                    ncbi_gene_id = record.ncbi_gene_id
                    mane_transcript_id = record.mane_transcript_id
                assert ncbi_gene_id == record.ncbi_gene_id
                assert mane_transcript_id == record.mane_transcript_id

            # NOTE(SW): HGNC IDs are only only available in GTF gene features
            hgnc_id = genes[ncbi_gene_id].hgnc_id

            for record in transcript_records:
                name = ';'.join((
                    record.symbol,
                    hgnc_id,
                    record.ncbi_gene_id,
                    record.mane_transcript_id,
                ))

                data = [
                    record.chrom,
                    record.start,
                    record.end,
                    name,
                    record.strand,
                ]
                print(*data, sep='\t', file=fh)


if __name__ == '__main__':
    main()
