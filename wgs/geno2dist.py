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
    p.add_argument('-o', '--outfile', default='outfile.dxy.txt')
    p.add_argument('infile')

    return p


def main( args ):

    geno2dxy( args )


def geno2dxy( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)
    lineparser.parse_grouping()

    cout('Grouping:')
    groups = lineparser.groups
    for k in lineparser.groups:
        cout(' %12s %3d' % (k, len(lineparser.groups[k])))

    group_keys = sorted(lineparser.groups.keys())
    cout(group_keys)

    # read whole genotype, and release all unused memory
    cerr('I: reading genotype file')
    haplotypes = lineparser.parse_np_haplotypes()

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



