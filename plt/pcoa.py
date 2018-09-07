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

    fig = plt.figure( figsize = (9, 9), dpi = args.dpi )

    fig_idx = 1
    colour_list = group_parser.colour_list()
    for pcx, pcy in combinations([0,1,2], 2):

        ax = fig.add_subplot(3, 1, fig_idx)
        fig_idx += 1

        make_plot(ax, pcoa[0][:,pcx], pcoa[0][:,pcy], colour_list)

    fig.tight_layout()
    fig.savefig(args.outfile)

    return

    import IPython
    IPython.embed()



    regions = df[df.columns[0]]
    region_boundaries = []
    start_idx, region_name = 0, regions[0]
    colours = cycle(colour_list[:3])
    for idx, region in enumerate(regions[1:]):
        if region != region_name:
            region_boundaries.append( (start_idx, idx-1, region_name, next(colours)) )
            start_idx, region_name = idx, region
    region_boundaries.append( (start_idx, idx-1, region_name, next(colours)) )

    print(region_boundaries)

    fig = plt.figure( figsize=(21, no_of_figures), dpi = args.dpi )

    for idx, c in enumerate(columns):
        cerr('I: plotting %s' % c)
        points = df[c]
        ax = fig.add_subplot(no_of_figures, 1, idx+1)
        make_plot(ax, points, region_boundaries, c)

    fig.tight_layout()
    fig.savefig(args.outfile)


def make_plot(axis, x, y, colours):
    axis.scatter( x, y, 0.2, c=colours )

def xxx_make_plot(axis, points, boundaries, label):
    for (start_idx, end_idx, region_name, region_colour) in boundaries:
        axis.scatter( np.arange(start_idx, end_idx), points[start_idx:end_idx], 0.25, c=region_colour)
        axis.set_ylim(0, 1.1)
        axis.set_xlim(-10, len(points) + 10)
        axis.get_xaxis().set_visible(False)
        axis.set_ylabel( label, fontsize=6 )






def haha():




    heights = df[df.columns[args.column - 1]]

    cerr('I: plotting...')
    #plt.bar( np.arange(0, len(heights)), heights, 1.0)
    #plt.plot( heights )
    plt.scatter( np.arange(0, len(heights)), heights, 0.05 )
    plt.savefig(args.outfile, dpi = args.dpi)

