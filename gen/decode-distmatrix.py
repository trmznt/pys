#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('decode (and optionally consolidate) encoded genetic distance')

    p.add_argument('-o', '--outfile', default='outdist.tsv',
                   help='output filename')
    p.add_argument('infiles', nargs='+',
                   help='input distance matrix files')

    return p


def decode_distmatrix(args):

    import pandas as pd
    import numpy as np
    from seqpy.core.bioio import tabutils
    from seqpy.core.funcs import distance

    cerr(f'[Reading from {len(args.infile)} file(s)]')
    dfs = [pd.read_table(infile) for infile in args.infiles]

    # do sanity checks
    columns = dfs[0].columns

    for df in dfs[1:]:
        if not columns.identical(df.columns):
            cexit('ERROR: headers do not match!')

    D = np.zeros(dfs[0].shape, dtype=np.int32)
    N = np.zeros(dfs[0].shape, dtype=np.int32)

    for df in dfs:
        d, n = distance.decode_distmatrix(df.values)
        D += d
        N += n

    distm = D/N
    distm_df = pd.DataFrame(data=distm, columns=columns)
    tabutils.write_file(args.outfile, distm_df)
    cerr(f'[Distance is written to {args.outfile}]')


def main(args):
    decode_distmatrix(args)


# EOF
