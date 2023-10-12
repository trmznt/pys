#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('join dataframes using a specific column as reference')

    p.add_argument('-o', '--outfile', default='outdata.tsv',
                   help='output filename')
    p.add_argument('infiles', nargs='+',
                   help='input dataframe files')

    return p


def join_dataframe(args):

    from seqpy.core.bioio import tabutils

    df = tabutils.read_file(args.infiles[0])
    for infile in args.infiles[1:]:

        if ':' in infile:
            infile, on_column = infile.split(':', 1)
        else:
            on_column = None

        right_df = tabutils.read_file(infile)
        if not on_column:
            on_column = right_df.columns[0]

        df = df.join(right_df.set_index(on_column), on=on_column, how='left')

    df.to_csv(args.outfile, sep='\t', index=False)
    cerr(f'Combined dataframe is written to {args.outfile}')


def main(args):
    join_dataframe(args)


# EOF
