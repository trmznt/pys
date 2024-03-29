#!/usr/bin/env spcli

import argparse
import numpy as np

from seqpy.core.sgk import sgio


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-b', '--bedfile', default='')
    p.add_argument('-s', '--samplefile', default='')
    p.add_argument('-o', '--outfile', default='outfile.txt')
    p.add_argument('-t', '--outtarget', default='')
    p.add_argument('--useGT', default=False, action='store_true',
                   help='Use GT field for genotype')
    p.add_argument('--mindepth', type=int, default=5)
    p.add_argument('--hetratio', type=float, default=-1,
                   help='The ratio of allele depth over total depth to call hets. '
                   'Value of 0.67 means if depth_of_major_allele/total_depth is < 0.67, '
                   'the genotype will be N.')
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

    if target_positions is None:
        # use all positions in source
        indexes = np.array(range(len(source_positions)))
        return indexes, source_positions

    for p in target_positions:
        try:
            indexes.append(source_positions.index((p[0], p[1])))
            real_targets.append(p)
        except ValueError:
            print(f'Warning: missing position: {p[0]} {p[1]}')
            indexes.append(-1)
            real_targets.append(p)

    return indexes, real_targets


def vcf2realmccoil(args):

    bed = None
    if args.bedfile:
        bed = read_bedfile(args.bedfile)
        if len(bed['regions']) > 0:
            print('Warning: the BED file contains region definition, but this tool only uses SNP definition')

    sample_set = None
    if args.samplefile:
        with open(args.samplefile) as f_sample:
            sample_set = set(filter(None, [x.strip() for x in f_sample]))
        print(f'Reading {len(sample_set)} samples from {args.samplefile}')

    ds = sgio.load_dataset(args.infile)

    samples = ds.sample_id.to_series()

    # SNPs
    snp_positions = list(zip(np.array(ds.contigs)[ds.variant_contig.values], ds.variant_position.values))
    if bed:
        snp_indexes, target_positions = prepare_snp_indexes(bed['snps'], snp_positions)
    else:
        snp_indexes, target_positions = prepare_snp_indexes(None, snp_positions)

    print(f'Reading {len(samples)} samples and {len(snp_positions)} SNPs from {args.infile}')

    if args.useGT:
        calldata_gt = ds.call_genotype
    else:
        calldata_ad = ds.call_AD
        _, d1, d2 = calldata_ad.shape
        allele_depths = calldata_ad.to_numpy()

    all_alleles = []
    selected_samples = []

    for n, s in enumerate(samples):

        if sample_set is not None and s not in sample_set:
            continue

        if args.useGT:

            # create alleles based on GT field
            alleles = calldata_gt[snp_indexes, n].sum(axis=1)
            alleles = alleles / 2
            alleles[alleles < 0] = -1

            all_alleles.append(alleles)
            selected_samples.append(s)
            continue

        sample_depths = allele_depths[snp_indexes, n]

        # add a priori 0.1 to total depths to avoid divide by zero (nan)
        total_depths = sample_depths.sum(axis=1, where=(sample_depths > 0)) + 0.1

        # get indexes for the highest reads per allele
        alleles = np.argmax(sample_depths, axis=1)

        # get indexes where total_depths are less than mindepth
        failed_masks = total_depths < args.mindepth

        # check heterozygosity
        if args.hetratio > 0:
            max_values = np.take_along_axis(
                sample_depths, np.expand_dims(alleles, axis=1), axis=1
            ).squeeze(axis=1)
            hetratios = max_values / total_depths
            het_masks = hetratios < args.hetratio
            alleles = alleles.astype(float)
            alleles[het_masks] = 0.5

        # change the positions where total depths are less than mindepths to X
        alleles[failed_masks] = -1

        # create the barcode string and add it barcodes container
        all_alleles.append(alleles)
        selected_samples.append(s)

        # import IPython; IPython.embed()
        # input()

    # write to text file
    with open(args.outfile, 'w') as fout:
        fout.write('ind\t' + '\t'.join(f'{p[0]}:{p[1]}' for p in target_positions) + '\n')
        for sample_code, alleles in zip(selected_samples, all_alleles):
            fout.write(f'{sample_code}\t')
            np.savetxt(fout, alleles, fmt='%.1f', delimiter='\t', newline='\t')
            fout.write('\n')

    print(f'TheRealMcCOIL (L={len(snp_indexes)}; N={len(selected_samples)}) written to {args.outfile}')

    # write target position
    if args.outtarget:
        with open(args.outtarget, 'w') as fout:
            for p in target_positions:
                fout.write(f'{p[0]}\t{p[1]}\t{p[2]}\n')
        print(f'Target positions written to {args.outtarget}')


def main(args):
    vcf2realmccoil(args)


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
