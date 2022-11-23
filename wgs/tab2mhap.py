#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile')
    p.add_argument('--codefile', required=True)
    p.add_argument('infile')
    return p


def tab2mhap(args):

    from seqpy.core.bioio.tabutils import read_file, write_file
    from seqpy.core.bioio.mhaputils import genotype_to_mhap

    geno_df = read_file(args.infile)
    mhcode_df = read_file(args.codefile)
    mhap_df, errs = genotype_to_mhap(geno_df, mhcode_df['MHCODE'])
    write_file(args.outfile, mhap_df)
    cerr(f'[microhaplotypes (N: {len(mhap_df)}, L: {len(mhap_df.columns) - 1}) '
         f'written to {args.outfile}]')
    if any(errs):
        cerr('Errors in generating microhaplotypes for:')
        cerr('\n'.join(errs))


def main(args):
    tab2mhap(args)

# EOF
