#!/usr/bin/env spcli

import argparse
import os

import pandas as pd

from seqpy import cerr
from seqpy.core.bioio import posutils


def init_argparser():
    p = argparse.ArgumentParser()
    p = posutils.init_argparser(p)

    p.add_argument('-o', '--outfile')
    p.add_argument('--outtarget', default='')
    p.add_argument('--samplefile', default='',
                   help='A headerless text file containing sample code per single line')
    p.add_argument('-m', '--metafile', default='',
                   help='Metafile containing columns to be added to output file, use eg. '
                   'metafile.tsv:SAMPLE,COUNTRY,REGION to add COUNTRY and REGION to output '
                   'using SAMPLE as index')
    p.add_argument('--useGT', default=False, action='store_true')
    p.add_argument('--mindepth', default=5, type=int,
                   help='Cut-off depth to be called missing variant, eg. mindepth = 5 '
                   'indicates variants with total depth < 5 will be mark as missing.')
    p.add_argument('--hetratio', default=0.67, type=float,
                   help='The ratio of allele depth over total depth to call hets. '
                   'Value of 0.67 means if depth_of_major_allele/total_depth is < 0.67, '
                   'the genotype will be N. Valid values are between 0.5 to 1.0. '
                   'Set to -1 for obtaining major allele')
    p.add_argument('--minaltdepth', default=2, type=int,
                   help='Threshold value for minor depth of a variant to be called heterozygote, '
                   'eg. minaltdepth = 2 indicates that variants with alternate reads >= 2 will be '
                   'marked as heterozygote, depending on the hetratio. Use hetratio = 0.999 if '
                   'hetratio is to be ignored.')
    p.add_argument('infile')

    return p


def zarr2tabular(args):

    from seqpy.core.sgk import sgio, sgutils
    from seqpy.core.bioio import tabutils

    # if providede, read posfile
    posdf = None
    if args.posfile:
        posdf = posutils.read_posfile(args=args)

    # load dataset
    ds = sgio.load_dataset(args.infile)

    # select SNPs
    if posdf is not None:
        ds = posdf.pos.sel_dataset(ds)

    # if need to select samples, performed here
    if args.samplefile:
        orig_N = ds.dims['samples']
        sample_df = pd.read_table(args.samplefile, header=None)
        ds = ds.sel(samples=ds.sample_id.isin(sample_df.iloc[:, 0].to_list()))
        curr_N = ds.dims['samples']
        if curr_N != len(sample_df):
            curr_samples = set(ds.sample_id.values)
            sel_samples = set(sample_df.iloc[:, 0])
            diff_samples = sel_samples - curr_samples
            raise ValueError(f'Samples not found: {diff_samples}')
        cerr(f'[Subsetting the samples from {orig_N} to {curr_N}]')

    # prepare metadata to be added to output file
    if args.metafile:
        samples = ds.sample_id.values
        sample_df, errs = tabutils.join_metafile(samples, args.metafile, percenttag=True)
    else:
        sample_df = None

    # convert using hetratio

    #import IPython; IPython.embed()

    cerr('[Converting alleles to tabular format...]')
    variants = sgutils.get_alleles(
        sgutils._allele_for_barcode,
        ds,
        hetratio=args.hetratio,
        mindepth=args.mindepth,
        useGT=args.useGT
    )

    tabular_df = tabutils.dataframe_from_variants(ds, variants)

    # saving to ouftile

    N, L = tabular_df.shape

    # add additional metadata, which already sorted based on samples
    if sample_df is not None:

        # assume 1st columns of tabular_df is SAMPLE, wchich can be replaced by sample_df
        tabular_df = pd.concat([sample_df, tabular_df.iloc[:, 1:]], axis=1)

    match ext := os.path.splitext(args.outfile)[1].lower():

        case '.tsv' | '.csv':
            tabular_df.to_csv(args.outfile, sep='\t' if ext == '.tsv' else ',', index=False)

        case '.txt':
            with open(args.outfile, 'w') as fout:
                fout.write('#BAR\tSAMPLE\n')
                for idx, row in tabular_df.iterrows():
                    fout.write(''.join(row[1:]) + '\t' + row[0] + '\n')

        case _:
            raise ValueError('other extension is not implemented yet')

    cerr(f'[Barcode (L={L-1};N={N}) is written to {args.outfile}]')

    if args.outtarget:
        target_posdf = posutils.posframe_from_dataset(ds)
        target_posdf.to_csv(args.outtarget, sep='\t', index=False)
        cerr(f'[Target position is written to {args.outtarget}]')


def main(args):
    zarr2tabular(args)

# EOF
