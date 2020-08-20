#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser


try:
    from matplotlib import pyplot as plt
except:
    cexit('ERR: require properly installed matplotlib')

import numpy as np

def init_argparser():
    p = arg_parser("Create depth plot")
    p.add_argument('-t', '--title', default='depth')
    p.add_argument('-d', '--mindepth', type=int, default=3)
    p.add_argument('--outgap', default='')
    p.add_argument('--statfile', default='')
    p.add_argument('-o', '--outfile', default='outplot.png')

    p.add_argument('infile')

    return p

def main( args ):

    depthplot( args )


import itertools
import yaml

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]

def depthplot( args ):

    # read data
    depth_list = []
    with open(args.infile) as fin:
        for line in fin:
            tokens = line.split()
            depth_list.append( (int(tokens[1]), int(tokens[2])) )

    # plot data
    length = depth_list[-1][0] + 1
    x = np.arange(length)
    y = np.zeros(length)

    for (pos, depth) in depth_list:
        y[pos] = depth

    # check for holes (depth < 1):
    holes = []
    for i in range(length):
        if y[i] < args.mindepth:
            holes.append(i)

    cerr('Missing region(s):')
    regions = []
    for region in list(ranges(holes)):
        regions.append(list(region))
        cerr(region)
    if args.outgap:
        yaml.dump({args.title: regions}, open(args.outgap, 'w'))

    if args.statfile:
        fig, axs = plt.subplots(1, 2, figsize=(30,5),
                gridspec_kw={'width_ratios': [5, 1]})
        ax = axs[0]
    else:
        fig, ax = plt.subplots(figsize=(25,5))

    ax.fill_between(x, y, facecolor="lightgreen", color="darkgreen", alpha=0.5)
    ax.set_yscale("log", base=10)
    ax.set_title(args.title)

    if args.statfile:

        # read data
        insert_size = []
        insert_count = []
        with open(args.statfile) as fin:
            for line in fin:
                if not line.startswith('IS'):
                    continue
                tokens = line.split()
                insert_size.append( int(tokens[1]) )
                insert_count.append( int(tokens[3]) )

        # plot data
        axs[1].hist(insert_size, len(insert_size), weights=insert_count)
        axs[1].set_yscale("log", base=10)

    fig.tight_layout()
    fig.savefig(args.outfile)


