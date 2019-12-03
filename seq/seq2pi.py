#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, grpparser, multisequence

import numpy as np
import itertools

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
    seqs = load(args.infile)

    # open and read group/meta file using groupfile
    cerr('[I - reading group information file]')
    group_parser = grpparser.GroupParser( args )
    group_parser.parse()

    group_seqs = {}

    for seq in seqs:
        grp = group_parser.group_info[seq.label.decode('ASCII')]
        if grp in group_seqs:
            group_seqs[grp].append( seq )
        else:
            ms = multisequence()
            ms.append( seq )
            group_seqs[grp] = ms

    print('Groups:')


def calc_pi( mseqs ):

    pi_array = []
    for (i,j) in itertools.combinations(range(len(mseqs)), 2):
        pi_array.append( calc_seq_diversity(mseqs[i], mseqs[j]))

    pi_array = np.array(pi_array)

    # return avg diversity and std dev

    return (,)


