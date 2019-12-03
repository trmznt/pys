#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser

import pandas


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile', default='outfile.txt')
    p.add_argument('-c', '--column', default=None)
    p.add_argument('infiles', nargs='+')
    return p


def main( args ):

    concat( args )


def concat( args ):

    columns = args.column.split(',') if args.column else None

    dfs = []
    for infile in args.infiles:
        df = pandas.read_table(infile)
        if columns:
            dfs.append( df.loc[:, columns ])
        else:
            dfs.append( df )

    # combine dfs
    new_df = pandas.concat( dfs, ignore_index=True)
    new_df.to_csv(args.outfile, '\t', index=False)






