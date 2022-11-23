#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile')
    p.add_argument('-s', '--samplefile')
    p.add_argument('infile')
    return p


def dist2dist(args):

    from seqpy.core.bioio.tabutils import read_file, write_file

    # read distance matrix file
    dist_df = read_file(args.infile)

    # read sample file
    sample_df = read_file(args.samplefile)

    samples = sample_df['SAMPLE']
    sample_idx = dist_df.columns.isin(samples)

    # sanity check
    if sum(sample_idx) != len(samples):
        raise ValueError(f'[Sample filtered: {sample_idx.sum()} is not identical '
                         f'with sample list: {len(samples)}')

    filtered_dist_df = dist_df.loc[sample_idx, sample_idx]
    write_file(args.outfile, filtered_dist_df)
    cerr(f'[New distance matrix is written to {args.outfile}]')


def main(args):
    dist2dist(args)

# EOF
