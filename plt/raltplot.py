#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def init_argparser(p=None):
    if not p:
        p = arg_parser("raltplot - plot ratio of alternate allele")
    p.add_argument('-o', '--outplot', default="outplot.ralt.png")
    p.add_argument('infile')

    return p

def main( args ):

    raltplot( args )


def raltplot( args ):

    infile = gzopen( args.infile )
    next(infile)        # read & discard header
    raltm = np.loadtxt(infile, delimiter='\t')

    ratios = []
    for r in raltm.flat:
        if 0.0 < r < 1.0:
            if r > 0.5:
                r = 1.0 - r
            ratios.append(r)

    ax = sns.distplot(ratios, bins=np.arange(0.01, 0.5, 0.01))
    plt.savefig(args.outplot)


    