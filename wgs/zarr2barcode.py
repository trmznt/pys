#!/usr/bin/env spcli

import argparse
import os

from seqpy import cerr
from seqpy.core.sgk import sgio, sgutils
from seqpy.core.bioio import posutils, tabutils


def init_argparser():
    p = argparse.ArgumentParser()
    p = posutils.init_argparser(p)

    p.add_argument('-o', '--outfile')
    p.add_argument('--outtarget', default='')
    p.add_argument('--useGT', default=False, action='store_true')
    p.add_argument('--mindepth', default=5, type=int)
    p.add_argument('--hetratio', default=0.67, type=float,
                   help='The ratio of allele depth over total depth to call hets. '
                   'Value of 0.67 means if depth_of_major_allele/total_depth is < 0.67, '
                   'the genotype will be N. Valid values are between 0.5 to 1.0. '
                   'Set to -1 for obtaining major allele')
    p.add_argument('infile')

    return p


def zarr2barcode(args):

    # if providede, read posfile
    posdf = None
    if args.posfile:
        posdf = posutils.read_posfile(args=args)

    # load dataset
    ds = sgio.load_dataset(args.infile)

    # select SNPs
    if posdf:
        ds = posdf.pos.sel_dataset(ds)

    # if need to select samples, performed here

    # convert using hetratio

    cerr('[Converting to barcode...]')
    variants = sgutils.get_alleles(
        sgutils._allele_for_barcode,
        ds,
        hetratio=args.hetratio,
        mindepth=args.mindepth,
        useGT=args.useGT
    )

    df_barcode = tabutils.dataframe_from_variants(ds, variants)

    # saving to ouftile
    match ext := os.path.splitext(args.outfile)[1].lower():

        case '.tsv' | '.csv':
            df_barcode.to_csv(args.outfile, sep='\t' if ext == '.tsv' else ',', index=False)

        case '.txt':
            with open(args.outfile, 'w') as fout:
                fout.write('#BAR\tSAMPLE\n')
                for idx, row in df_barcode.iterrows():
                    fout.write(''.join(row[1:]) + '\t' + row[0] + '\n')

        case _:
            raise ValueError('other extension is not implemented yet')

    N, L = df_barcode.shape
    cerr(f'[Barcode (L={L-1};N={N}) is written to {args.outfile}]')


def main(args):
    zarr2barcode(args)

# EOF
