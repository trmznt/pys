#!/usr/bin/env spcli


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser

from itertools import cycle, combinations
from matplotlib import pyplot as plt

try:
    from matplotlib import pyplot as plt
    import allel
except:
    cexit('ERR: require properly installed matplotlib and scikit-allel')
import numpy as np

def init_argparser():
    p = arg_parser("Create PCoA plot based on distance matrix file")
    p = grpparser.init_argparser(p)
    p.add_argument('--dpi', type=int, default=600)
    p.add_argument('--dotsize', type=float, default=0.25)
    p.add_argument('-o', '--outfile', default="outplot.pcoa.png")

    p.add_argument('infile')

    return p


def main( args ):

    pcoa( args )


def pcoa( args ):

    cerr('I: reading group info')
    group_parser = grpparser.GroupParser( args )
    group_parser.parse()

    with open(args.infile, 'rb') as infile:

        cerr('I: reading sample header...')
        samples = next(infile).decode('UTF-8').strip().split()
        groups = group_parser.assign_groups(samples)

        cerr('I: reading distance matrix')
        distm = np.loadtxt(infile, delimiter='\t')

    pcoa = allel.pcoa(distm)

    fig = plt.figure( figsize = (27, 9), dpi = args.dpi )

    fig_idx = 1
    colour_list = group_parser.colour_list()
    for pcx, pcy in combinations([0,1,2], 2):

        ax = fig.add_subplot(1, 3, fig_idx)
        fig_idx += 1

        make_plot(ax, pcoa[0][:,pcx], pcoa[0][:,pcy], colour_list, args.dotsize)

    fig.tight_layout()
    fig.savefig(args.outfile)


def make_plot(axis, x, y, colours, dotsize):
    axis.scatter( x, y, dotsize, c=colours )

