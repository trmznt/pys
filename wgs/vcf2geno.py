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
            fields=['samples', 'variants/CHROM', 'variants/POS', 
                'calldata/AD'])

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
    cerr('[I: writing position file]')
    with open(pos_file, 'w') as outfile:
        outfile.write('CHROM\tPOS\n')
        for (chrom, pos) in zip(
                vcfset['variants/CHROM'], vcfset['variants/POS']):
            print(chrom, pos)
            outfile.write(chrom)
            outfile.write('\t')
            outfile.write(str(pos))
            outfile.write('\n')

    # write genotype by converting the genotype
    cerr('[I: writing genotype file]')
    with open(geno_file, 'w') as outfile:
        outfile.write( '\t'.join( list( vcfset['samples'])))
        outfile.write('\n')
        for gts in vcfset['calldata/AD']:
            #np.savetxt( outfile, majgeno(gts), fmt="%d", delimiter="\t" )
            outfile.write( '\t'.join( '%d' % x for x in majgeno(gts)) )
            outfile.write('\n')


def majgeno(genotypes):

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

def majgeno_misshet(genotypes):

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


def hetgeno(genotypes, threshold):

    pass

