#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import fastdx

import numpy as np

def init_argparser(p=None):

    p = naltparser.init_argparser( with_group=False )
    p.add_argument('--includepos', default='')
    p.add_argument('--revfmt', default=False, action='store_true' )
    p.add_argument('-o', '--outfile', default='outfile.haplotypes.txt')

    return p


def main( args ):

    nalt2haplotypes( args )


def nalt2haplotypes( args ):

    nalt_parser = naltparser.NAltLineParser( args, with_group=False )
    region = nalt_parser.parse_whole()
    samples = nalt_parser.parse_samples()

    if args.includepos:
        with open(args.includepos) as f:
            poslines = [ x.strip().split() for x in f ]
            if poslines[0][0] == 'CHROM' and poslines[0][1] == 'POS':
                del poslines[0]
    else:
        poslines = None

    if poslines:
        region.filter_poslines(poslines, inplace=True)

    # read whole genotype, and release all unused memory
    cerr('[I - converting to haplotypes]')
    haplotypes = region.haplotypes()
    positions = region.P

    lines = []
    for i in range(len(haplotypes)):
        snps = []
        for j in range(len(haplotypes[i])):
            allel = haplotypes[i][j]
            if allel == 0:
                snps.append( positions[j][2])
            elif allel == 2:
                snps.append( positions[j][3])
            elif allel < 0:
                snps.append( 'X')
            else:
                snps.append( 'N' )
        lines.append( (samples[i], ''.join(snps)))

    if args.outfile:
        with open(args.outfile, 'w') as f:
            if args.revfmt:
                for line in lines:
                    f.write( '{}\t{}\n'.format(line[1], line[0]))
            else:
                for line in lines:
                    f.write( '{}\t{}\n'.format(line[0], line[1]))

