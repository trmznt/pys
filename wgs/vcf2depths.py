#!/usr/bin/env python3

import argparse
import allel
import numpy as np
import pandas as pd


# this script generate a dataframe file for allele depths for each samples

def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outfile', default='outdepths.txt')
    p.add_argument('-b', '--bedfile', default='')
    p.add_argument('--minlen', type=int, default=10,
                   help='Cutoff for minimum region length.')
    p.add_argument('--maxinterval', type=int, default=5,
                   help='Cutoff for the interval still considered as a single spanning region.')
    p.add_argument('infile')
    return p


def read_bedfile(bedfile):
    """ return a dictionary containing SNPs positions and region positions """
    bed = {'snps': [], 'regions': []}
    with open(bedfile) as f:
        for line in f:
            tokens = line.strip().split()
            start_pos, end_pos = int(tokens[1]), int(tokens[2])
            span = end_pos - start_pos
            if span <= 0:
                raise ValueError(f'BED file improper value: {tokens[0]}, {start_pos}, {end_pos}')
            if end_pos - start_pos == 1:
                bed['snps'].append((tokens[0], end_pos, tokens[3]))
            elif end_pos - start_pos > 1:
                bed['regions'].append((tokens[0], start_pos, end_pos, tokens[3]))

    return bed


def create_aggregate_from_bed(positions, bedinfo):
    aggregate_positions = []
    positions = list(positions)
    for bed_item in bedinfo['snps']:
        try:
            i = positions.index((bed_item[0], bed_item[1]))
            aggregate_positions.append((bed_item[0], bed_item[1], 1, i, i + 1))
        except ValueError:
            print(f'Warning: missing position {bed_item}')

    for bed_item in bedinfo['regions']:
        # import IPython; IPython.embed()
        # input()
        chr, start_pos, end_pos = bed_item[0], bed_item[1] + 1, bed_item[2]
        try:
            start_i = positions.index((chr, start_pos))
            end_i = positions.index((chr, end_pos))

            # sanity check
            length = end_pos - start_pos
            if positions[start_i + length][1] != end_pos:
                print(f'Warning: inconsistent region for {bed_item}: missing {length-(end_i-start_i)} positions')

            aggregate_positions.append((chr, start_pos, length, start_i, end_i))

            # import IPython; IPython.embed()

        except ValueError:
            print(f'Warning: missing region {bed_item}')

    return aggregate_positions


def create_aggregate_positions(positions, delta_pos=1):
    """ return [(chrom, position, length, index_start, index_end), ...] """

    print(f'Warning: maxinterval = {delta_pos}')
    aggregate_positions = []
    curr_c = curr_p = start_p = start_i = None
    for i, (c, p) in enumerate(positions):
        if not curr_c:
            curr_c = c
            curr_p = start_p = p
            start_i = i
            continue
        if c == curr_c and p < curr_p:
            raise ValueError('SNP are not sorted by positions')
        if (p - curr_p > delta_pos) or curr_c != c:
            aggregate_positions.append((curr_c, start_p, curr_p - start_p + 1, start_i, i))
            curr_c = c
            curr_p = start_p = p
            start_i = i

        curr_p = p

    aggregate_positions.append((curr_c, start_p, p - start_p + 1, start_i, i))

    return aggregate_positions


def vcf2depths(args):

    bed = None
    if args.bedfile:
        bed = read_bedfile(args.bedfile)

    vcf = allel.read_vcf(args.infile, fields=['samples', 'variants', 'calldata/DP', 'calldata/GT', 'calldata/AD'])

    if bed:
        aggregate_positions = create_aggregate_from_bed(zip(vcf['variants/CHROM'], vcf['variants/POS']), bed)
    else:
        aggregate_positions = create_aggregate_positions(zip(vcf['variants/CHROM'], vcf['variants/POS']),
                                                         delta_pos=args.maxinterval)

        # filter for length > minlen
        aggregate_positions = [p for p in aggregate_positions if p[2] >= args.minlen]

    samples = vcf['samples']

    allele_depths = allel.GenotypeArray(vcf['calldata/AD'])
    L, N = len(aggregate_positions), len(samples)

    # region_depths has row for position and column for samples
    region_depths = np.zeros((L, N), int)

    for n_idx in range(N):
        print(f'\rsample: {n_idx}', end='')
        for l_idx in range(L):
            pos_info = aggregate_positions[l_idx]
            region = allele_depths[pos_info[3]:pos_info[4], n_idx]
            depths = region.sum(axis=1, where=(region > 0))
            if pos_info[2] == 1:
                region_depths[l_idx, n_idx] = depths.min()
            else:
                region_depths[l_idx, n_idx] = np.quantile(depths[~vcf['variants/INDEL'][pos_info[3]:pos_info[4]]], 0.25)
    print()

    df = pd.DataFrame(region_depths, columns=samples, index=[f'{p[0]}:{p[1]}+{p[2]}' for p in aggregate_positions])
    df.to_csv(args.outfile, sep='\t')


def main(args):
    vcf2depths(args)


# this script can be run independently or through seqpy spcli
if __name__ == '__main__':
    args = init_argparser().parse_args()
    if args.debug:
        from ipdb import launch_ipdb_on_exception
        with launch_ipdb_on_exception:
            main(args)
    else:
        main(args)

# EOF
