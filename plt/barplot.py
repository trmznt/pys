#!/usr/bin/env spcli


from seqpy import cout, cerr
from seqpy.cmds import arg_parser


try:
    from matplotlib import pyplot as plt
    import pandas
except:
    cexit('ERR: require properly installed matplotlib and pandas')
import numpy as np

def init_argparser():
    p = arg_parser("Create bar plot from columnar data")
    p.add_argument('--column', type=int, required=True)
    p.add_argument('--asc', action="store_true")
    p.add_argument('--desc', action="store_true")
    p.add_argument('--xlabel', default='')
    p.add_argument('--ylabel', default='')
    p.add_argument('--title', default='')
    p.add_argument('--dpi', type=int, default=600)
    p.add_argument('-o', '--outfile', default="outplot.png")

    p.add_argument('infile')

    return p


def main( args ):

    barplot( args )


def barplot( args ):

    cerr('I: reading data...')
    df = pandas.read_table( args.infile )

    column = df.columns[args.column - 1]
    cerr('I: selecting column %s' % column)

    if args.asc:
        cerr('I: sorting ascending...')
        df = df.sort_values( column )
    elif args.desc:
        cerr('I: sorting descending...')
        df = df.sort_values( column, ascending=False)

    heights = df[column]

    cerr('I: plotting...')
    #plt.bar( np.arange(0, len(heights)), heights, 1.0)
    #plt.plot( heights )
    plt.scatter( np.arange(0, len(heights)), heights, 0.25 )
    if args.xlabel:
        plt.xlabel(args.xlabel)
    if args.ylabel:
        plt.ylabel(args.ylabel)
    if args.title:
        plt.title(args.title)
    plt.savefig(args.outfile, dpi = args.dpi)

