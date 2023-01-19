#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile')
    p.add_argument('--codefile', required=True)
    p.add_argument('-m', '--metafile', default='',
                   help='Metafile containing columns to be added to output file, use eg. '
                   'metafile.tsv:SAMPLE,COUNTRY,REGION to add COUNTRY and REGION to output '
                   'using SAMPLE as index')
    p.add_argument('infile')
    return p


def tab2mhap(args):

    import pandas as pd
    from seqpy.core.bioio.tabutils import read_file, write_file, join_metafile
    from seqpy.core.bioio.mhaputils import genotype_to_mhap

    geno_df = read_file(args.infile)
    mhcode_df = read_file(args.codefile)
    mhap_df, errs = genotype_to_mhap(geno_df, mhcode_df['MHCODE'])

    N, L = len(mhap_df), len(mhap_df.columns) - 1

    # prepare metadata to be added to output file
    if args.metafile:
        samples = mhap_df['SAMPLE']
        sample_df, meta_errs = join_metafile(samples, args.metafile, percenttag=True)
    else:
        sample_df = None

    # add additional metadata, which already sorted based on samples
    if sample_df is not None:

        # assume 1st columns of tabular_df is SAMPLE, wchich can be replaced by sample_df
        mhap_df = pd.concat([sample_df, mhap_df.iloc[:, 1:]], axis=1)

    write_file(args.outfile, mhap_df)
    cerr(f'[microhaplotypes (N: {N}, L: {L}) '
         f'written to {args.outfile}]')
    if any(errs):
        cerr('Errors in generating microhaplotypes for:')
        cerr('\n'.join(errs))


def main(args):
    tab2mhap(args)

# EOF
