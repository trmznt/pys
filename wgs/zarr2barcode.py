#!/usr/bin/env spcli

import argparse
import numpy as np

from seqpy.core.sgk import sgio, sgutils
from seqpy.core.bioio import posutils


def init_argparser():
    p = argparse.ArgumentParser()
    p = posutils.init_argparser(p)

    p.add_argument('-o', '--outfile')
    p.add_argument('--outtarget', default='')
    p.add_argument('--useGT', default=False, action='store_true')
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

    variants = []
    for var_idx in ds.variants:

        variants.append(
            sgutils.get_alleles(
                sgutils._allele_for_barcode,
                ds,
                hetratio=args.hetratio,
                mindepth=args.mindepth,
                useGT=args.useGT
            )
        )

    # saving to outfile


def main(args):
    zarr2barcode(args)

# EOF
