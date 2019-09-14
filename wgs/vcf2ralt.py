#!/usr/bin/env spcli

# this script convert VCF file to 1) ralt (ratio of alternate allele) matrix and 2) SNP/position manifest
# and 3) nmdp (no of minor read). Please note that this script needs to be run under spcli (seqpy)

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
import time

try:
    import allel
except:
    cerr('ERR: require properly installed scikit-allel!')
import numpy as np

from seqpy.core.cfuncs import genoutils

def init_argparser():
    p = arg_parser("Convert VCF to ratio of alternate ref dataset")
    p.add_argument('--autofilename', default=False, action='store_true')
    p.add_argument('--mindepth', default=1, type=int)
    p.add_argument('-o', '--outfile', default='outdata')
    p.add_argument('infile')

    return p

def main(args):
    vcf2ralt( args )


def vcf2ralt( args ):

    # read the mighty VCF file
    cerr('[I: reading VCF...]')
    start_time = time.monotonic()
    vcfset = allel.read_vcf(args.infile,
            fields=['samples', 'variants/CHROM', 'variants/POS', 'variants/REF',
                'variants/ALT', 'variants/SNPEFF_GENE_NAME',
                'variants/SNPEFF_AMINO_ACID_CHANGE', 'calldata/AD'])
    cerr('[I: read %s site, %s samples in %d secs]' % (len(vcfset['variants/CHROM']),
         len(vcfset['samples']), time.monotonic() - start_time))

    if args.autofilename:
        args.outfile = 'r-%d-%d' % (
                len(vcfset['samples']), len(vcfset['variants/POS'])
        )

    pos_file = args.outfile + '.pos.txt'
    geno_file = args.outfile + '.ralt.txt'
    nmdp_file = args.outfile + '.nmdp.txt'

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

    # write genotype by converting the genotype
    cerr('[I: writing genotype file]')
    with open(geno_file, 'w') as outfile, open(nmdp_file, 'w') as outnmdp:
        outfile.write( '\t'.join( list( vcfset['samples'])))
        outfile.write('\n')
        outnmdp.write( '\t'.join( list( vcfset['samples'])))
        outnmdp.write('\n')
        c = 0
        for gts in vcfset['calldata/AD']:
            #np.savetxt( outfile, majgeno(gts), fmt="%d", delimiter="\t" )
            ratios, nread = genoutils.ralt(gts)
            #outfile.write( '\t'.join( '%4.3f' % x for x in genoutils.ralt(gts)) )
            #outfile.write('\n')
            np.savetxt(outfile, [ratios], delimiter='\t', fmt='%4.3f')
            np.savetxt(outnmdp, [nread], delimiter='\t', fmt='%d')
            c += 1
            if c % 1000 == 0:
                cerr('[I: writing site %d]' % c)

    cerr('[I: finish writing %d sites in %d secs]' % (c, time.monotonic() - start_time))


def ralt(genotypes, mindepth = 1):

    data = np.zeros( shape = len(genotypes), dtype=float )
    for idx, gt in enumerate(genotypes):
        gt_tot = gt[0] + gt[1]
        if gt_tot < mindepth:
            data[idx] = -1.0
            continue
        data[idx] = gt[0]/gt_tot

    return data
