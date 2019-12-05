#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, grpparser, multisequence
from seqpy.core.funcs import profiles

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

    seq2fst( args )


def seq2fst( args ):

    # open and read sequence file
    cerr('[I - reading sequence file %s]' % args.infile)
    seqs = load(args.infile)

    # open and read group/meta file using groupfile/metafile if available
    if args.groupfile or args.metafile:
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
    else:
        cexit('[ERR - seq2fst.py requires group information!]')

    print('Groups:')
    for g in group_seqs:
        avg, stddev = calc_pi(group_seqs[g])
        print('  %20s [%3d]: %f %f' % (g, len(group_seqs[g]), avg, stddev))


def calc_fst( mseqs1, mseqs2 ):

    pi_array = []
    for (i,j) in itertools.combinations(range(len(mseqs)), 2):
        pi_array.append( calc_propdist(mseqs[i], mseqs[j]))

    pi_array = np.array(pi_array)

    # return avg diversity and std dev
    return (np.mean(pi_array), np.std(pi_array))


def to_genotype_array(*mseqs):

    # genotype array is array of site vs sample vs [0,0]
    # 0 for the major allele, 1 for the minor allele

    # create matrix profile first
    full_mseqs = multisequence()
    for mseq in mseqs:
    	full_mseqs.extend( mseq )

    na_profiles = profile.na_profiles( full_mseqs )
    



def calc_propdist(seq1, seq2):

    if len(seq1) != len(seq2):
        raise RuntimeError('[ERR: seq %s and %s do not have similar length' % (seq1, seq2))

    p = l = 0.0
    for i in range(len(seq1)):
        if seq1[i] == b'X' or seq2[i] == b'X':
            continue
        if seq1[i] == seq2[i]:
            l += 1
        elif seq1[i] == b'N' or seq2[i] == b'N':
            p += 0.5
            l += 1
        else:
            p += 1
            l += 1

    return p/l




