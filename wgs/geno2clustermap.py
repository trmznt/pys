#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create a clustermap from the genotype file


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-c', default=None)
    p.add_argument('--outplot', default='clustremap.png')

    return p


def main( args ):

    geno2clustermap( args )


def geno2clustermap( args ):

    genoparser = tabparser.GenotypeLineParser( args )
    np_haplotypes = genoparser.parse_np_haplotypes()

    with open(args.c) as colourfile:
    	next(colourfile)
    	colors = [ x.strip().split()[1] for x in colourfile ]
    axis = sb.clustermap(np_haplotypes, row_colors = colors, yticklabels=False, method='average')
    plt.show()


