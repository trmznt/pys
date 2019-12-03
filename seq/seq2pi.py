#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, grpparser

import numpy as np

def init_argparser(p=None):

    if p is None:
        p = arg_parser('seq2pi.py - calculate nucleotide diversity (pi)')

    p = grpparser.init_argparser( p )
    p.add_argument('-o', '--outfile', default='outfile.pi.txt')
    p.add_argument('infile')

    return p


def main( args ):

    seq2pi( args )


def seq2pi( args ):

    # open and read sequence file
    cerr('[I - reading sequence file %s]' % args.infile)
    seq = load(args.infile)

    # open and read group/meta file using groupfile
