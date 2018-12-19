#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import io


def init_argparser(p=None):
    if not p:
        p = arg_parser("raltplot - plot ratio of alternate allele")
    p.add_argument('-o', '--outplot', default="outplot.ralt.png")
    p.add_argument('ralt_file')
    p.add_argument('nmdp_file')

    return p

def main( args ):

    raltplot( args )


def raltplot( args ):

    cerr('[I - reading input files]')

    ratios = []
    n_mdps = []

    with gzopen( args.ralt_file ) as ralt_file, gzopen( args.nmdp_file) as nmdp_file:
        next(ralt_file)        # read & discard header
        next(nmdp_file)
        #raltm = np.loadtxt(ralt_file, delimiter='\t')
        #nmdpm = np.loadtxt(nmdp_file, delimiter='\t')
        c = 0
        for ralt_line, nmdp_line in zip(ralt_file, nmdp_file):
            raltm = np.loadtxt(io.StringIO(ralt_line), dtype=np.float, delimiter='\t')
            nmdpm = np.loadtxt(io.StringIO(nmdp_line), dtype=np.short, delimiter='\t')

    #cerr('[I - analyzing input files]')

            for i in range(len(raltm)):
                r = raltm[i]
                n = nmdpm[i]
                if 0.0 < r < 1.0:
                    if r > 0.5:
                        r = 1.0 - r
                    ratios.append(r)
                    n_mdps.append(n)

            c += 1
            if c % 500 == 0:
                cerr('[I - reading site %d]' % c)

    cerr('[I - finish reading %d site]' % c)
    cerr('[I - plotting %d het SNPs]' % len(ratios))
    g = sns.jointplot(ratios, n_mdps, s=0.5).set_axis_labels("ratio", "minor DP")
    #ax = sns.jointplot(x='r', y='n', data=df)
    g.ax_joint.set_yscale('log')

    #ax = sns.distplot(ratios, bins=np.arange(0.01, 0.5, 0.01))
    plt.savefig(args.outplot)
