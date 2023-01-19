#!/usr/bin/env spcli

import argparse
from seqpy import cerr


def init_argparser():

    p = argparse.ArgumentParser()

    p.add_argument('-o', '--outfile')
    p.add_argument('--components', type=int, default=3,
                   help='Number of principal components to be written, default is 3')
    p.add_argument('infile')
    return p


def tab2mca(args):

    from seqpy.core.bioio import tabutils
    import pandas as pd
    import prince
    import time

    cerr(f'Reading genomic data from {args.infile}')
    df = tabutils.read_file(args.infile)
    df.geno.set_alleles(missings=['X'])
    allele_df = df.geno.get_alleles()
    sample_df = df.geno.get_samples()
    metacolumns_df = df.geno.get_metacolumns()

    mca = prince.MCA(
        n_components=args.components,
        n_iter=6,
        copy=True,
        check_input=True,
        engine='auto',
        random_state=int(time.time())
    )

    N, L = allele_df.shape
    cerr(f'[Processing {N} samples with {L} markers]')
    mca = mca.fit(allele_df)
    CA = mca.row_coordinates(allele_df)
    CA.rename(columns={i: f'PC{i+1}' for i in range(args.components)},
              inplace=True)

    # prepare metadata to be added to output file
    joined_df = pd.concat([metacolumns_df, CA], axis=1)

    if args.outfile:
        tabutils.write_file(args.outfile, joined_df)
        cerr(f'[MCA written to {args.outfile} with {args.components} components]')


def main(args):
    tab2mca(args)

# EOF
