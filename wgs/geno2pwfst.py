#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import itertools
import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.pwfst.txt')
    p.add_argument('infile')

    return p


def main( args ):

    geno2pwfst( args )


def geno2pwfst( args ):
    """ perform pair-wise FST by population """

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.diploid_translator)
    lineparser.parse_grouping()

    cout('Grouping:')
    groups = lineparser.groups
    for k in lineparser.groups:
        cout(' %12s %3d' % (k, len(lineparser.groups[k])))

    FST = [] # FST indexed by group_keys
    group_keys = sorted(lineparser.groups.keys())

    # read whole genotype, and release all unused memory
    cerr('I: reading genotype file')
    allel_array = lineparser.parse_all()
    cerr('I: generating genotype array')
    genoarray = allel.GenotypeArray( allel_array )
    del allel_array

    cerr('I: counting alleles')
    ac = {}
    for g in group_keys:
        ac[g] = genoarray.count_alleles( subpop = groups[g])

    cerr('I: calculating FST')
    M = np.zeros( (len(group_keys), len(group_keys)) )
    for (i,j) in itertools.permutations( range(len(group_keys)), 2):

        i_group = group_keys[i]
        j_group = group_keys[j]
        fst, _, _, _ = allel.stats.blockwise_hudson_fst(
                ac[i_group], ac[j_group], blen=10)
        M[i,j] = M[j,i] = fst

    with open(args.outfile, 'wt') as outfile:
        # write header:
        outfile.write('%s\n' % ('\t'.join(group_keys)))
        np.savetxt(outfile, M, delimiter='\t')


    return
