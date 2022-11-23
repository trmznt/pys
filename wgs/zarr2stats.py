#!/usr/bin/env spcli

import argparse

import sgkit as sg

from seqpy import cerr
from seqpy.core.sgk import sgio, sgutils, sgstats
from seqpy.core.bioio import posutils, tabutils


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('--biallelic', action='store_true', default=False)
    p.add_argument('--biallelic_ratio_threshold', type=float, default=-1)
    p.add_argument('--indel_ratio_cutoff', type=float, default=-1)
    p.add_argument('infile')

    return p


def stat_biallele(ds, biallelic_ratio_threshold=-1, indel_ratio_cutoff=1):

    linebuf = []
    _ = linebuf.append

    ds = sgstats.variant_stats(ds)

    # count number of strict biallelic SNPs
    cerr('[Calculating strict biallelic variants]')
    strict_biallelic = (ds['variant_biallelic_ratio'] >= 1.0)
    _(f'Number of strict biallelic variant: {strict_biallelic.sum().values}')

    cerr('[Calculating strict SNV]')
    strict_SNV = (ds['variant_indel_ratio'] <= 0.0)
    _(f'Number of strict single nucleotide variant: {strict_SNV.sum().values}')

    cerr('[Calculating strict biallelic SNV]')
    strict_biallelic_SNV = (strict_biallelic & strict_SNV)
    _(f'Number of strict biallelic SNV: {strict_biallelic_SNV.sum().values}')

    # count number of relax biallelic SNPs
    if biallelic_ratio_threshold > 0:
        cerr('[Calculating relax biallelic variants]')
        relax_biallelic = (ds['variant_biallelic_ratio'] >= biallelic_ratio_threshold).sum().values
        _(f'Number of relax biallelic variant: {relax_biallelic}')

    if indel_ratio_cutoff > 0:
        cerr('[Calculating relax SNV]')
        relax_SNV = (ds['variant_indel_ratio'] <= indel_ratio_cutoff).sum().values
        _(f'Number of relax single nucleotide variant: {relax_SNV}')

    if biallelic_ratio_threshold > 0 and indel_ratio_cutoff > 0:
        cerr('[Calculating relax biallelic SNV]')
        relax_biallelic_SNV = (
            (ds['variant_biallelic_ratio'] >= biallelic_ratio_threshold) &
            (ds['variant_indel_ratio'] <= indel_ratio_cutoff)
        ).sum().values
        _(f'Number of relax biallelic SNV: {relax_biallelic_SNV}')

    # print out text
    cerr('\n'.join(linebuf))


def zarr2stats(args):

    ds = sgio.load_dataset(args.infile, 7)

    if args.biallelic:
        stat_biallele(ds, args.biallelic_ratio_threshold, args.indel_ratio_cutoff)


def main(args):
    zarr2stats(args)

# EOF
