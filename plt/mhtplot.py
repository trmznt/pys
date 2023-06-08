#!/usr/bin/env spcli


from seqpy import cout, cerr
from seqpy.cmds import arg_parser

from itertools import cycle
from matplotlib import pyplot as plt

try:
    from matplotlib import pyplot as plt
    import pandas
except:
    cexit('ERR: require properly installed matplotlib and pandas')
import numpy as np

def init_argparser():
    p = arg_parser("Create multiple manhattan plot from columnar data")
    p.add_argument('--column', default=None)
    p.add_argument('--dpi', type=int, default=600)
    p.add_argument('--dotsize', type=float, default=0.25)
    p.add_argument('--autoyscale', default=False, action='store_true')
    p.add_argument('-o', '--outfile', default="outplot.png")

    p.add_argument('infile')

    return p


def main( args ):

    mhtplot( args )

colour_list = [ '#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928',
                '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']


def mhtplot( args ):

    cerr('I: reading data...')
    df = pandas.read_table( args.infile )

    if args.column:
        columns = [ df.columns[int(i)] for i in args.column.split(',') ]
    else:
        columns = df.columns[3:]

    no_of_figures = len(columns)

    regions = df[df.columns[0]]
    region_boundaries = []
    start_idx, region_name = 0, regions[0]
    colours = cycle(colour_list[:5])
    for idx, region in enumerate(regions):
        if region != region_name:
            region_boundaries.append( (start_idx, idx-1, region_name, next(colours)) )
            start_idx, region_name = idx, region
    region_boundaries.append( (start_idx, idx, region_name, next(colours)) )

    print(region_boundaries)

    fig = plt.figure( figsize=(21, no_of_figures * 5), dpi = args.dpi )

    for idx, c in enumerate(columns):
        cerr('I: plotting %s' % c)
        points = df[c]
        ax = fig.add_subplot(no_of_figures, 1, idx+1)
        make_plot(ax, points, region_boundaries, c, args.dotsize, autoyscale=args.autoyscale)

    fig.tight_layout()
    fig.savefig(args.outfile)


def make_plot(axis, points, boundaries, label, dotsize=0.25, autoyscale=False):
    for (start_idx, end_idx, region_name, region_colour) in boundaries:
        end_idx += 1
        axis.scatter( np.arange(start_idx, end_idx), points[start_idx:end_idx], dotsize, c=region_colour)
        if not autoyscale:
            axis.set_ylim(0, 1.1)
        axis.set_xlim(-10, len(points) + 10)
        axis.get_xaxis().set_visible(False)
        axis.set_ylabel( label, fontsize=6 )
