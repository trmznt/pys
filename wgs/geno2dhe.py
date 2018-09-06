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

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.dhe.txt')
    p.add_argument('infile')

    return p


def main( args ):

    geno2dhe( args )


def geno2dhe( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)
    lineparser.parse_grouping()

    cout('Grouping:')
    groups = lineparser.groups
    for k in lineparser.groups:
        cout(' %12s %3d' % (k, len(lineparser.groups[k])))

    FST = [] # FST indexed by group_keys
    group_keys = sorted(lineparser.groups.keys())
    cout(group_keys)

    # read whole genotype, and release all unused memory
    cerr('I: reading genotype file')
    allel_array = lineparser.parse_all()
    cerr('I: generating genotype array')
    genoarray = allel.GenotypeArray( allel_array )
    del allel_array

    cerr('I: calculating He')
    He = 1 - np.sum( genoarray.count_alleles().to_frequencies()**2, axis=1 )

    He_groups = {}
    pHe = None
    for g in groups:
        He_groups[g] = 1 - np.sum(
                genoarray.count_alleles(subpop=groups[g]).to_frequencies()**2,
                axis = 1)
        if pHe is None:
            pHe = He_groups[g]
        else:
            pHe = pHe + He_groups[g]

    dHe = He - pHe/len(He_groups)

    cerr('I: writing output file')
    with open(args.outfile, 'wt') as outfile:
        outfile.write('CHROM\tPOS\tdHe\tHe\t%s\n' % '\t'.join(group_keys))

        for i in range(len(He)):
            posinfo = lineparser.position[i]

            outfile.write('%s\t%s\t%5.4f\t%5.4f\t%s\n' % (
                posinfo[0], posinfo[1], dHe[i], He[i], 
                '\t'.join( '%5.4f' % He_groups[g][i] for g in group_keys))
            )

