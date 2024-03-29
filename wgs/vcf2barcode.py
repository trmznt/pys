#!/usr/bin/env python3

import argparse
import os
import allel
import numpy as np


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-b', '--bedfile', required=True)
    p.add_argument('-o', '--outfile', default='outbarcode.txt')
    p.add_argument('-t', '--outtarget', default='')
    p.add_argument('--mindepth', type=int, default=5)
    p.add_argument('--hetratio', type=float, default=-1,
                   help='The ratio of allele depth over total depth to call hets. '
                   'Value of 0.67 means if depth_of_major_allele/total_depth is < 0.67, '
                   'the genotype will be N.')
    p.add_argument('--usevcfpos', default=False, action='store_true',
                   help='Use vcf for outtarget, ie chrom, pos, ref, alt')
    p.add_argument('--debug', default=False, action='store_true')
    p.add_argument('infile')
    return p


def read_bedfile(bedfile):
    """ return a dictionary containing SNPs positions and region positions """
    bed = {'snps': [], 'regions': []}
    with open(bedfile) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tokens = line.strip().split()
            if len(tokens) == 2:
                bed['snps'].append((tokens[0], int(tokens[1]), 'snp'))
                continue
            start_pos, end_pos = int(tokens[1]), int(tokens[2])
            span = end_pos - start_pos
            if span <= 0:
                raise ValueError(f'BED file improper value: {tokens[0]}, {start_pos}, {end_pos}')
            if span == 1:
                bed['snps'].append((tokens[0], end_pos, tokens[3]))
            elif span > 1:
                bed['regions'].append((tokens[0], start_pos, end_pos, tokens[3]))

    return bed


def prepare_snp_indexes(target_positions, source_positions, snp_refalt):
    indexes = []
    real_targets = []

    # check duplication in source_positions:
    if len(source_positions) != len(set(source_positions)):
        raise ValueError('Duplicate position found, please normalize or filter the vcf file first')

    for p in target_positions:
        try:
            idx = source_positions.index((p[0], p[1]))
            indexes.append(idx)
            if snp_refalt is not None:
                real_targets.append((source_positions[idx][0], source_positions[idx][1],
                                    snp_refalt[idx][0], snp_refalt[idx][1][0]))
            else:
                real_targets.append(p)
        except ValueError:
            print(f'Warning: missing position: {p[0]} {p[1]}')
            indexes.append(-1)
            if snp_refalt is not None:
                real_targets.append((p[0], p[1], 'X', 'X'))
            real_targets.append(p)

    return indexes, real_targets


def vcf2barcode(args):

    bed = read_bedfile(args.bedfile)
    if len(bed['regions']) > 0:
        print('Warning: the BED file contains region definition, but this tool only uses SNP definition')

    vcf = allel.read_vcf(args.infile, fields=['samples', 'variants', 'calldata/DP', 'calldata/GT', 'calldata/AD'])

    samples = vcf['samples']

    # SNPs
    snp_positions = list(zip(vcf['variants/CHROM'], vcf['variants/POS']))
    snp_refalt = list(zip(vcf['variants/REF'], vcf['variants/ALT'])) if args.usevcfpos else None
    snp_indexes, target_positions = prepare_snp_indexes(bed['snps'], snp_positions, snp_refalt)

    # adding mmissing data at the end of all array data to accomdate missing positions, ie. -1 indexing position
    calldata_ad = vcf['calldata/AD']
    _, d1, d2 = calldata_ad.shape
    allele_depths = allel.GenotypeArray(np.append(calldata_ad, np.full((1, d1, d2), -1, dtype=np.int16), axis=0))
    refs = np.append(vcf['variants/REF'], ['X'])[snp_indexes]
    alts = np.append(vcf['variants/ALT'], [['X', '', '']], axis=0)[snp_indexes]
    alleles = np.insert(alts, 0, refs, axis=1)
    N, L = len(samples), len(snp_indexes)
    barcodes = []

    for n in range(N):
        sample_depths = allele_depths[snp_indexes, n]

        # add a priori 0.1 to total depths to avoid divide by zero (nan)
        total_depths = sample_depths.sum(axis=1, where=(sample_depths > 0)) + 0.1

        # get indexes for the highest reads per allele
        indexes_max = np.argmax(sample_depths, axis=1)

        # get indexes where total_depths are less than mindepth
        failed_masks = total_depths < args.mindepth

        # get alleles based on indexes_max
        bar = np.take_along_axis(alleles, np.expand_dims(indexes_max, axis=1), axis=1).squeeze(axis=1)

        # check heterozygosity
        if args.hetratio > 0:
            max_values = np.take_along_axis(
                sample_depths, np.expand_dims(indexes_max, axis=1), axis=1
            ).squeeze(axis=1)
            hetratios = max_values / total_depths
            het_masks = hetratios < args.hetratio
            bar[het_masks] = 'N'

        # change the positions where total depths are less than mindepths to X
        bar[failed_masks] = 'X'

        # create the barcode string and add it barcodes container
        barcodes.append(''.join(bar))

        # import IPython; IPython.embed()
        # input()

    # write to text file
    match os.path.splitext(args.outfile)[1]:

        case '.txt':
            with open(args.outfile, 'w') as fout:
                for sample_code, barcode in zip(samples, barcodes):
                    fout.write(f'{sample_code}\t{barcode}\n')

        case '.csv' | '.tsv' as ext:
            sep = ',' if ext == '.csv' else '\t'
            header = ['SAMPLES'] + [f'{p[0]}:{p[1]}' for p in target_positions]
            with open(args.outfile, 'w') as fout:
                # write header
                fout.write(sep.join(header))
                fout.write('\n')
                # write content
                for sample_code, barcode in zip(samples, barcodes):
                    fout.write(sample_code)
                    fout.write(sep)
                    fout.write(sep.join(barcodes))
                    fout.write('\n')

        case _ as ext:
            raise RuntimeError(f'Please use extention .txt, .csv or .tsv in output filename, instead of {ext}')

    print(f'Barcodes (L={len(snp_indexes)}) written to {args.outfile}')

    # write target position
    if args.outtarget:
        with open(args.outtarget, 'w') as fout:
            if args.usevcfpos:
                fout.write('CHROM\tPOS\tREF\tALT\n')
            for p in target_positions:
                if args.usevcfpos:
                    fout.write(f'{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\n')
                else:
                    fout.write(f'{p[0]}\t{p[1]}\t{p[2]}\n')
        print(f'Target positions written to {args.outtarget}')


def main(args):
    vcf2barcode(args)


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
