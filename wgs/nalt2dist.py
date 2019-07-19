#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import fastdx

import numpy as np

def init_argparser(p=None):

    p = naltparser.init_argparser()
    p.add_argument('--countmissing', default=False, action='store_true')
    p.add_argument('-o', '--outfile', default='outfile.dist.txt')

    return p


def main( args ):

    nalt2dist( args )


def nalt2dist( args ):

    nalt_parser = naltparser.NAltLineParser( args, with_group=False, with_position=False )
    whole_region = nalt_parser.parse_whole()


    # read whole genotype, and release all unused memory
    cerr('[I - converting to haplotypes]')
    haplotypes = whole_region.haplotypes()

    cerr('[I - calculating pairwise dxy for %d samples]' % len(haplotypes))
    distm = fastdx.pwdistm( haplotypes, args.countmissing )

    cerr('[I - writing to %s]' % args.outfile)
    with open(args.outfile, 'w') as outfile:
        outfile.write( '\t'.join(nalt_parser.parse_samples()) )
        outfile.write( '\n')

        # write the matrix
        np.savetxt(outfile, distm, delimiter='\t', fmt='%.6f')
