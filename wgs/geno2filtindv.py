#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# remove samples with certain missingness


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np
import itertools

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('--cutoff', type=float, default=0.0)
    p.add_argument('--outfile', default='outfile.geno.txt')

    return p


def main( args ):

    geno2filtindv( args )


def geno2filtindv( args ):

    cerr('I: reading genotype file')
    genoparser = tabparser.GenotypeLineParser( args )
    cerr('I: generating haplotypes')
    haplotypes = genoparser.parse_haplotypes()

    cerr('I: scanning haplotypes')
    flags = [True] * len(haplotypes)
    for idx, haplo in enumerate(haplotypes):
        missingness = haplo.count(b'-')/len(haplo)
        if missingness > args.cutoff:
            flags[idx] = False
        #cerr('I: %4d : %f' % (idx, haplo.count(b'-')/len(haplo))

    cerr('I: filtering samples')
    outfile = open(args.outfile, 'w')
    genoparser2 = tabparser.GenotypeLineParser( args )
    samples = itertools.compress(genoparser2.samples, flags)
    outfile.write('\t'.join(samples)); outfile.write('\n')

    for posline, genoline in genoparser2.parse_raw_lines():
        new_genotypes = itertools.compress( genoline.strip().split(), flags)
        outfile.write('\t'.join(new_genotypes)); outfile.write('\n')

    cerr('I: writing for %d samples' % flags.count(True))
    