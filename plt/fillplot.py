#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('fillplot')
    p.add_argument('-o', '--outfile', default='outplot.png')
    p.add_argument('-x', required=True,
                   help="Column name whose values to be used for x-axis")
    p.add_argument('-y', required=True,
                   help="Comma separated column names whose values to be used for y-axis")
    p.add_argument('infile')

    return p


def fillplot(args):

    try:
        from matplotlib import pyplot as plt

    except ImportError:
        cexit('ERR: require properly installed matplotlib')

    from seqpy.core.bioio import tabutils
    import itertools

    df = tabutils.read_file(args.infile)

    fig, ax = plt.subplots(figsize=(25, 5))
    colors = itertools.cycle(['forestgreen', 'royalblue', 'lightcoral'])

    x = df[args.x]
    for column in args.y.split(','):

        y = df[column]
        color = next(colors)
        ax.fill_between(x, y, facecolor=color, color=color, alpha=0.5)

    fig.tight_layout()
    fig.savefig(args.outfile)
    cerr(f'Plot is written to {args.outfile}')


def main(args):

    fillplot(args)

# EOF
