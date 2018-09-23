#!/usr/bin/env spcli

# this script convert VCF file to three file:
# 1) genotype matrix, 2) sample manifest and 3) SNP/position manifest
# please note that this script needs to be run under spcli (seqpy)

from seqpy import cout, cerr
from seqpy.cmds import arg_parser



try:
    import allel
except:
    cerr('ERR: require properly installed scikit-allel!')
import numpy as np

def init_argparser():
    p = arg_parser("Convert VCF to genotype dataset")
    p.add_argument('-o', '--outfile', default='outdata')
    p.add_argument('--filter', default='')
    p.add_argument('--majority', default=False, action='store_true')
    p.add_argument('--minhetratio', default=0.25, type=float)
    p.add_argument('infile')

    return p

def main(args):
    vcf2geno( args )


def vcf2geno( args ):

    # read the mighty VCF file
    cerr('[I: reading VCF...]')
    vcfset = allel.read_vcf(args.infile,
            fields=['samples', 'variants/CHROM', 'variants/POS', 'variants/REF',
                'variants/ALT', 'variants/SNPEFF_GENE_NAME',
                'variants/SNPEFF_AMINO_ACID_CHANGE', 'calldata/AD'])
    cerr('[I: read %s site, %s samples]' % (len(vcfset['variants/CHROM']),
         len(vcfset['samples'])))

    sample_file = args.outfile + '.indv.txt'
    pos_file = args.outfile + '.pos.txt'
    geno_file = args.outfile + '.geno.txt'

    # write samples
    cerr('[I: writing sample file]')
    with open(sample_file, 'w') as outfile:
        outfile.write('SAMPLE\n')
        for indv_code in vcfset['samples']:
            outfile.write(indv_code)
            outfile.write('\n')

    # write position
    # position file is:
    # CHROM POS REF ALT
    cerr('[I: writing position file]')
    with open(pos_file, 'w') as outfile:
        outfile.write('CHROM\tPOS\tREF\tALT\tGENE\tAACHANGE\n')
        for (chrom, pos, ref, alt, gene, aachange) in zip(
                vcfset['variants/CHROM'], vcfset['variants/POS'],
                vcfset['variants/REF'], vcfset['variants/ALT'],
                vcfset['variants/SNPEFF_GENE_NAME'],
                vcfset['variants/SNPEFF_AMINO_ACID_CHANGE']):
            outfile.write('%s\t%d\t%s\t%s\t%s\t%s\n' %
                (chrom, pos, ref, alt[0], gene, aachange))

    # set scheme
    if args.majority:
        func = majgeno
        threshold, r_threshold = 0, 0 
    else:
        func = hetgeno
        threshold, r_threshold = args.minhetratio, 1.0 - args.minhetratio

    # write genotype by converting the genotype
    cerr('[I: writing genotype file]')
    with open(geno_file, 'w') as outfile:
        outfile.write( '\t'.join( list( vcfset['samples'])))
        outfile.write('\n')
        for gts in vcfset['calldata/AD']:
            #np.savetxt( outfile, majgeno(gts), fmt="%d", delimiter="\t" )
            outfile.write( '\t'.join( '%d' % x for x in func(gts, threshold, r_threshold)) )
            outfile.write('\n')


def majgeno(genotypes, threshold, r_threshold):

    data = np.zeros( shape = len(genotypes), dtype=int )
    for idx, gt in enumerate(genotypes):
        gt_tot = gt[0] + gt[1]
        if gt_tot == 0:
            data[idx] = -1
        elif gt[0] > gt[1]:
            data[idx] = 0
        else:
            data[idx] = 2

    return data


def majgeno_misshet(genotypes, threshold, r_threshold):

    data = np.zeros( shape = len(genotypes), dtype=int )
    for idx, gt in enumerate(genotypes):
        gt_tot = gt[0] + gt[1]
        if gt_tot == 0:
            data[idx] = 1
        elif gt[0] > gt[1]:
            data[idx] = 0
        else:
            data[idx] = 2

    return data


def hetgeno(genotypes, threshold, r_threshold):

    data = np.zeros( shape = len(genotypes), dtype=int )
    for idx, gt in enumerate(genotypes):
        gt_tot = gt[0] + gt[1]
        if gt_tot == 0:
            data[idx] = -1
            continue
        ratio = gt[0]/gt_tot
        if ratio < threshold:
            data[idx] = 0
        elif ratio > r_threshold:
            data[idx] = 2
        else:
            data[idx] = 1

    return data
