#!/usr/bin/env python3
import argparse
import pathlib


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel_fp', required=True, type=pathlib.Path)
    parser.add_argument('--hartwig_panel_fp', required=True, type=pathlib.Path)

    args = parser.parse_args()

    if not args.panel_fp.exists():
        parser.error(f'Input file {args.panel_fp} does not exist')
    if not args.hartwig_panel_fp.exists():
        parser.error(f'Input file {args.hartwig_panel_fp} does not exist')

    return args


def main():
    # Get command line arguments
    args = get_arguments()

    # Read in data
    hartwig_panel = get_panel_data(args.hartwig_panel_fp, 'gene')
    umccr_panel = get_panel_data(args.panel_fp, 'ensembl_gene_symbol')

    # Prepare then format data
    umccr_panel = prepare_data(umccr_panel, hartwig_panel)
    umccr_formatted = format_data(umccr_panel)

    # Output data
    header_tokens = list(hartwig_panel.values())[0].keys()
    print(*header_tokens, sep='\t')

    for data in umccr_formatted:
        print(*data, sep='\t')


def prepare_data(umccr_panel, hartwig_panel):
    # For UMCCR genes with ambiguous role, use the role set by Hartwig where available otherwise
    # set role as tsgene for now. Some UMCCR genes will have no role, also set role as tsgene for now
    #
    # NOTE(SW): using tsgene as default is fine for now since it only impacts the likelihood
    # calculation, which we are not currently using in our curation process
    #
    # Records are mutated inplace
    for record in umccr_panel.values():
        role_none = record['oncogene'] in {'NA', 'FALSE'} and record['tsgene'] in {'NA', 'FALSE'}
        role_ambiguous = record['oncogene'] == 'TRUE' and record['tsgene'] == 'TRUE'
        assert not (role_none and role_ambiguous)

        # Skip resolved gene role
        if not (role_none or role_ambiguous):
            continue

        # Default to tsgene if no further information available from Hartwig
        if record['ensembl_gene_symbol'] not in hartwig_panel:
            record['oncogene'] = 'FALSE'
            record['tsgene'] = 'TRUE'
            continue

        # Set both TRUE here and target one FALSE below
        record['oncogene'] = 'TRUE'
        record['tsgene'] = 'TRUE'

        # Apply information from Hartwig
        hartwig_gene_role = hartwig_panel[record['ensembl_gene_symbol']]['likelihoodType']
        if hartwig_gene_role == 'ONCO':
            record['tsgene'] = 'FALSE'
        elif hartwig_gene_role == 'TSG':
            record['oncogene'] = 'FALSE'
        else:
            assert False

    # NOTE(SW): returning for clarity
    return umccr_panel


def format_data(umccr_panel):
    reportable_somatic_columns = ['true'] * 7
    reportable_germline_columns = ['ANY'] * 4

    data = list()
    for record in umccr_panel.values():
        assert not (record['oncogene'] == 'TRUE' and record['tsgene'] == 'TRUE')
        if record['oncogene'] == 'TRUE':
            likelihood_type = 'ONCO'
        elif record['tsgene'] == 'TRUE':
            likelihood_type = 'TSG'
        else:
            assert False

        data.append((
            record['ensembl_gene_symbol'],
            *reportable_somatic_columns,
            likelihood_type,
            *reportable_germline_columns,
            '',
            'false',
        ))

    return data


def get_panel_data(fp, key_name):
    data = dict()
    with fp.open('r') as fh:
        header_tokens = fh.readline().rstrip().split('\t')
        for line in fh:
            record = {k: v for k, v in zip(header_tokens, line.rstrip().split('\t'))}

            # NOTE(SW): skip records that are not present in Ensembl 105
            # Currently: TRA, TRB, IGH, IGK, IGL
            if record[key_name] == 'NA':
                continue

            assert record[key_name] not in data
            data[record[key_name]] = record
    return data


if __name__ == '__main__':
    main()
