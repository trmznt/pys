#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# report heterozygosity statistics parameters


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outpath', default='outfile.hetstats')

    return p


def main( args ):

    geno2hetstats( args )


def geno2hetstats( args ):

    genoparser = tabparser.GenotypeLineParser( args )
    genoparser.set_translator(tabparser.GenotypeLineParser.haploid_translator)

    samples = genoparser.parse_sample()
    N = len(samples)
    sample_hets = np.zeros( N )
    sample_miss = np.zeros( N )
    n = np.arange( N )
    c = 0

    cerr('[I: reading genotype file...]')
    outsite = open( args.outpath + '.site-het.txt', 'w')
    outsite.write('CHROM\tPOS\tREGION\tN_HETS\tHETS\tN_MISS\tMISS\n')
    for posinfo, genovector in genoparser.parse():
        c += 1
        site_hets = 0
        site_miss = 0
        for i in n:
            g = genovector[i]
            if g == 1:
                sample_hets[i] += 1
                site_hets += 1
            elif g == -1:
                sample_miss[i] += 1
                site_miss += 1

        if c % 1000 == 0:
            cerr('[I: reading site %d - %s %s]' % (c, posinfo[0], posinfo[1]))
        outsite.write('%s\t%s\t%s\t%d\t%5.4f\t%d\t%5.4f\n' % 
            (    posinfo[0], posinfo[1], posinfo[4],
                site_hets, site_hets/(N-site_miss), site_miss, site_miss/N)
        )

    sample_typed = c - sample_miss
    het_ratios = sample_hets/sample_typed
    miss_ratios = sample_miss/c
    outindv = open( args.outpath + '.indv-het.txt', 'w')
    outindv.write('SAMPLE\tN_HETS\tHETS\tN_MISS\tMISS\n')
    for i in n:
        outindv.write('%s\t%d\t%5.4f\t%d\t%5.4f\n' % 
            (samples[i], sample_hets[i], het_ratios[i], sample_miss[i], miss_ratios[i])
        )


