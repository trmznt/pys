#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('zarr2dist - calculate proportional genetic distance')

    p.add_argument('--encode', default=False, action='store_true',
                   help="encode distance with Szudzik's pairing function")
    p.add_argument('--mindepth', default=5, type=int,
                   help='minimum depth per variant')
    p.add_argument('--allelenumber', default=2, type=int,
                   help='max number of allele (use 2 for biallelic) [2]')
    p.add_argument('-o', '--outfile', default='outdist.tsv',
                   help="file output [outdist.tsv]")
    p.add_argument('infile')

    return p


def zarr2dist(args):

    import pandas as pd
    from seqpy.core.bioio import tabutils
    from seqpy.core.sgk.sgio import load_dataset
    from seqpy.core.funcs import distance

    cerr(f'[Reading genotype data from {args.infile}]')
    ds = load_dataset(args.infile)

    cerr('[Generating major allele matrix...]')
    alleles = distance.AD_to_alleles(ds.call_AD.values, args.mindepth, args.allelenumber)

    cerr('[Calculating pairwise distances...]')
    distm = distance.pairwise_distances(alleles)

    if not args.encode:
        # split distm to proportional distance
        cerr('[Decoding distance matrix...]')
        d, N = distance.decode_distmatrix(distm)
        distm = d / N

    distm_df = pd.DataFrame(data=distm, columns=ds.sample_id)
    tabutils.write_file(args.outfile, distm_df)
    cerr(f'[Distance is written to {args.outfile}]')


def main(args):
    zarr2dist(args)

# EOF
