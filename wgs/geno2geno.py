#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import itertools
import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', required=True)

    return p


def main( args ):

    geno2geno( args )


def geno2geno( args ):
    """ perform pair-wise FST by population """

    genoparser = tabparser.GenotypeLineParser( args )

    sample_header = genoparser.get_sample_header()
    filename = args.outfile
    with open(filename+'.geno.txt','wt') as outfile, open(filename+'.pos.txt','wt') as outpos:

        outfile.write( sample_header )
        outfile.write('\n')
        outpos.write( genoparser.get_position_header() )
        outpos.write('\n')
        c = 0

        for posline, genoline in genoparser.parse_raw_lines():

            # split posline
            tokens = posline.split()
            if (tokens[0], tokens[1]) in genoparser.include_positions:
                outfile.write( genoline )
                outpos.write( posline )
                c += 1

    cerr('I: writing %d positions' % c)




