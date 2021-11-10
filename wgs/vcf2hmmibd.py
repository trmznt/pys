#!/usr/bin/env python3

import argparse
import allel
import numpy as np
import pandas as pd


# this script generate hmmIBD input file from a VCF file

def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outfile', default='hmmibd-genotype.txt')
    p.add_argument('-d', '--mindepth', type=int, default=1)
    p.add_argument('-t', '--translationfile', default='')
    p.add_argument('--debug', default=False, action='store_true')
    p.add_argument('infile')
    return p


def read_chrom_translation(infile):

    d = {}
    f = open(infile)
    for line in f:
        tokens = line.strip().split()
        d[tokens[0]] = int(tokens[1])
    return d

def vcf2hmmibd(args):

    print(f'Reading VCF file: {args.infile}')
    vcf = allel.read_vcf(args.infile, fields=['samples', 'variants', 'calldata/DP', 'calldata/GT', 'calldata/AD'])

    samples = vcf['samples']

    # convert chrom name to integer
    trans_d = read_chrom_translation(args.translationfile)
    chrom = [trans_d[x] for x in vcf['variants/CHROM']]
    coordinates = np.column_stack((chrom, vcf['variants/POS']))

    # create genotypes
    allele_depths = allel.GenotypeArray(vcf['calldata/AD'])
    genotypes = np.argmax(allele_depths, axis=2)

    # set for missing values or lower than mindepth
    total_depths = allele_depths.sum(axis=2, where=(allele_depths > 0))
    genotypes[total_depths < args.mindepth] = -1

    columns = ['chrom', 'pos'] + list(samples)
    genotypes = np.hstack((coordinates, genotypes))

    df = pd.DataFrame(genotypes, columns=columns)
    df.to_csv(args.outfile, sep='\t', index=False)
    print(f'Input file for hmmIBD written at: {args.outfile}')


def main(args):
    vcf2hmmibd(args)


# this script can be run independently or through seqpy spcli
if __name__ == '__main__':
    args = init_argparser().parse_args()
    if args.debug:
        from ipdb import launch_ipdb_on_exception
        with launch_ipdb_on_exception():
            main(args)
    else:
        main(args)

# EOF
