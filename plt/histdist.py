#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def init_argparser():
    p = arg_parser("Create PCoA plot based on distance matrix file")
    p.add_argument('--dpi', type=int, default=600)
    p.add_argument('-o', '--outplot', default="outplot.pcoa.png")

    p.add_argument('infile')

    return p


def main( args ):

    histdist( args )


def histdist( args ):

    with open(args.infile, 'rb') as infile:

        cerr('I: reading sample header...')
        samples = next(infile).decode('UTF-8').strip().split()
        
        cerr('I: reading distance matrix')
        distm = np.loadtxt(infile, delimiter='\t')

    d1 = distm.ravel()
    cerr('[I: total pairwising %d]' % len(d1))
    cerr('[I: max=%f]' % max(d1))
    ax = sns.distplot(distm.ravel())
    plt.show()




