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
    p.add_argument('-o', '--outpath', default='outfile')

    return p


def main( args ):

    geno2filtsite( args )


def geno2filtsite( args ):

    cerr('I: reading genotype file')
    genoparser = tabparser.GenotypeLineParser( args )

    outgeno = open(args.outpath + '.geno.txt', 'w')
    outgeno.write(genoparser.get_sample_header()); outgeno.write('\n')

    outpos = open(args.outpath + '.pos.txt', 'w')
    outpos.write(genoparser.get_position_header()); outpos.write('\n')

    c = 0
    for posline, genoline in genoparser.parse_raw_lines():
        tokens = genoline.split()
        missingness = tokens.count('-1') / len(tokens)
        if missingness > args.cutoff:
            continue
        outgeno.write(genoline)
        outpos.write(posline)
        c += 1

    cerr('I: writing for %d sites' % c)
