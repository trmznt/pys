#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import itertools
import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('--posindex', required=True)
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('infile')

    return p


def main( args ):

    geno2geno( args )


def geno2geno( args ):
    """ perform pair-wise FST by population """

    lineparser = tabparser.GenotypeLineParser( args )
    
