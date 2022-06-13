#!/usr/bin/env spcli

import argparse
import numpy as np
import pandas as pd

from seqpy.core.sgk import sgio


# this script generate hmmIBD input file from a VCF file
# use major genotype for hets call

def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outfile', default='hmmibd-genotype.txt')
    p.add_argument('-d', '--mindepth', type=int, default=1)
    p.add_argument('-s', '--samplefile', default='')
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

    sample_set = None
    if args.samplefile:
        with open(args.samplefile) as f_sample:
            sample_set = set(filter(None, [x.strip() for x in f_sample]))
        print(f'Reading {len(sample_set)} samples from {args.samplefile}')

    ds = sgio.load_dataset(args.infile)

    samples = ds.sample_id.to_series()
    collected_samples = []
    indices = []

    if sample_set:
        for i, s in samples.iteritems():
            if s in sample_set:
                collected_samples.append(s)
                indices.append(i)

    # convert chrom name to integer
    trans_d = read_chrom_translation(args.translationfile)
    chrom = [trans_d[x] for x in np.array(ds.contigs)[ds.variant_contig]]
    coordinates = np.column_stack((chrom, ds.variant_position))

    #import IPython; IPython.embed()

    # create genotypes
    allele_depths = ds.call_AD
    if sample_set:
        allele_depths = allele_depths[:, indices, :]

    genotypes = ds.call_AD.argmax(axis=2).compute()
    
    # set for missing values or lower than mindepth
    total_depths = allele_depths.to_numpy().sum(axis=2, where=(allele_depths > 0))
    genotypes.values[total_depths < args.mindepth] = -1

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
