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
import scipy.spatial

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.dxy.txt')
    p.add_argument('infile')

    return p


def main( args ):

    geno2dxy( args )


def geno2dxy( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)
    lineparser.parse_grouping()

    cout('Grouping:')
    groups = lineparser.groups
    for k in lineparser.groups:
        cout(' %12s %3d' % (k, len(lineparser.groups[k])))

    group_keys = sorted(lineparser.groups.keys())
    cout(group_keys)

    # read whole genotype, and release all unused memory
    cerr('I: reading genotype file')
    allel_array = lineparser.parse_all()
    cerr('I: generating genotype array')
    genoarray = allel.GenotypeArray( allel_array )
    del allel_array

    cerr('I: generating allele count array')
    gac = genoarray.to_allele_counts()

    cerr('I: calculating pairwise dxy')
    c_distm = allel.pairwise_dxy(range(len(gac)), gac)
    distm = scipy.spatial.distance.squareform(c_distm)

    #import IPython
    #IPython.embed()

    cerr('I: writing to outfile')
    with open(args.outfile, 'wb') as outfile:
        outfile.write( lineparser.get_sample_header(True) )
        outfile.write( b'\n')

        # write the matrix
        np.savetxt(outfile, distm, delimiter='\t') #, fmt='%.5f')
