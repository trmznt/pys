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
    p = arg_parser("Report het proportion per sample")
    p.add_argument('-o', '--outfile', default='outhet.txt')

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--minhetratio', default=0.25, type=float)
    g.add_argument('--usegt', default=False, action='store_true')

    p.add_argument('infile')

    return p

def main(args):
    vcf2hetstats( args )


def vcf2hetstats( args ):

    # read the mighty VCF file
    cerr('[I: reading VCF...]')
    vcfset = allel.read_vcf(args.infile,
            fields=['samples', 'variants/CHROM', 'variants/POS', 'calldata/AD', 'calldata/GT'])

    hets = []

    if args.usegt:
        GT = vcfset['calldata/GT']
        assert len(GT[0]) == len (vcfset['samples'])
        for i, sample in enumerate(vcfset['samples']):
            gts = GT[:,i]
            count = het = 0
            for gt in gts:
                if gt[0] == -1:
                    continue
                if gt[0] != gt[1]:
                    het += 1
                count += 1

            hets.append( (sample, count, het, het/count))

    else:
        cerr('ERR: can not proceed!')


    with open(args.outfile, 'w') as outfile:
        outfile.write('SAMPLE\tN_SNPS\tN_HETS\tF_HETS\n')
        for (sample, n_snp, n_het, f_het) in hets:
            outfile.write('%s\t%d\t%d\t%5.4f\n' % (sample, n_snp, n_het, f_het))

    cerr('Output written to %s' % args.outfile)
