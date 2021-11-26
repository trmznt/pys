#!/usr/bin/env python3

import argparse
import allel
import numpy as np


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-b', '--bedfile', required=True)
    p.add_argument('-o', '--outfile', default='outtable.txt')
    p.add_argument('-t', '--outtarget', default='')
    p.add_argument('-f', '--field', default='GT')
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
            start_pos, end_pos = int(tokens[1]), int(tokens[2])
            span = end_pos - start_pos
            if span <= 0:
                raise ValueError(f'BED file improper value: {tokens[0]}, {start_pos}, {end_pos}')
            if end_pos - start_pos == 1:
                bed['snps'].append((tokens[0], end_pos, tokens[3]))
            elif end_pos - start_pos > 1:
                bed['regions'].append((tokens[0], start_pos, end_pos, tokens[3]))

    return bed


def prepare_snp_indexes(target_positions, source_positions):
    indexes = []
    real_targets = []

    # check duplication in source_positions:
    if len(source_positions) != len(set(source_positions)):
        raise ValueError('Duplicate position found, please normalize or filter the vcf file first')

    for p in target_positions:
        try:
            indexes.append(source_positions.index((p[0], p[1])))
            real_targets.append(p)
        except ValueError:
            print(f'Warning: missing position: {p[0]} {p[1]}')
            indexes.append(-1)
            real_targets.append(p)

    return indexes, real_targets


def vcf2barcode(args):

    bed = read_bedfile(args.bedfile)
    if len(bed['regions']) > 0:
        print('Warning: the BED file contains region definition, but this tool only uses SNP definition')

    calldata_field = f'calldata/{args.field.upper()}'
    vcf = allel.read_vcf(args.infile, fields=['samples', 'variants', calldata_field])

    samples = vcf['samples']

    # SNPs
    snp_positions = list(zip(vcf['variants/CHROM'], vcf['variants/POS']))
    snp_indexes, target_positions = prepare_snp_indexes(bed['snps'], snp_positions)

    # get the relevant data
    calldata = vcf[calldata_field]
    # assert 0
    if args.field == 'GT':
        genotype_array = allel.GenotypeArray(calldata)
    else:
        genotype_array = calldata

    N = len(samples)
    chroms = vcf['variants/CHROM'][snp_indexes]
    positions = vcf['variants/POS'][snp_indexes]

    # straight write to file
    outf = open(args.outfile, 'w')
    outf.write('CHROM\t' + '\t'.join(chroms) + '\n')
    outf.write('POS\t' + '\t'.join(str(p) for p in positions) + '\n')

    for n in range(N):
        outf.write(samples[n] + '\t')
        if args.field == 'GT':
            ga = genotype_array[:,n][[snp_indexes]]
            _, items = ga.get_display_items(len(snp_indexes), len(snp_indexes))
            outf.write('\t'.join(items))
        else:
            np.savetxt(outf, genotype_array[:,n][snp_indexes], fmt='%d', newline='\t')
        outf.write('\n')

    outf.close()
    print(f'Table written to {args.outfile}')


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
