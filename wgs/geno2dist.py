#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an He file with the following format:
# dHe He POP1 POP2 POP3 ...


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np
from itertools import combinations

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.dist.txt')

    return p


def main( args ):

    geno2dist( args )


def geno2dist( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)

    # read whole genotype, and release all unused memory
    cerr('I: reading genotype file')
    haplotypes = lineparser.parse_haplotypes()

    cerr('I: calculating pairwise dxy')
    distm = pairwise_dxy( haplotypes )

    cerr('I: writing to outfile')
    with open(args.outfile, 'wb') as outfile:
        outfile.write( lineparser.get_sample_header(True) )
        outfile.write( b'\n')

        # write the matrix
        np.savetxt(outfile, distm, delimiter='\t', fmt='%.6f')



def pairwise_dxy(haplotypes):

    n = len(haplotypes)
    distm = np.zeros((n,n))
    idxs = list(range(len(haplotypes[0])))

    for i,j in combinations( range(n), 2):
        cerr('I: pairwising %d - %d' % (i,j))
        x = haplotypes[i]
        y = haplotypes[j]
        d = 0
        c = 0
        for idx in idxs:
            xi = x[idx]
            yi = y[idx]
            if xi == b'-' or yi == b'-': continue
            if xi != yi: d += 1
            c += 1
        distm[i,j] = distm[j,i] = d/c

    return distm



