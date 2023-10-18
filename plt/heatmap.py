#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('creating a heatmap from a Dataframe file')
    p.add_argument('-o', '--outplot')
    p.add_argument('-c', '--colorfile',
                   help=("a TSV/CSV file containg the group, COLOR and MARKER columns "
                         "for defining the color and marker of each data point"))
    p.add_argument('--hue', default=None,
                   help="Column name to used to differentiate data points")
    p.add_argument('--dpi', type=int, default=600)
    # p.add_argument('--jitter', type=float, default=0.005)
    p.add_argument('--equal_aspect', default=False, action='store_true',
                   help='set both axis to equal aspect (scale)')
    p.add_argument('--fontscale', type=float, default=1,
                   help="font scale, default=1")
    p.add_argument('-s', '--size', type=int, default=10,
                   help="size of markers")
    p.add_argument('-w', '--width', type=float, default=8.0,
                   help="Width & height of each scatter plot, in inch. Default is 8.")
    p.add_argument('--legend_out', action='store_true', default=False,
                   help="Legend will be put outside (right of) plot")
    p.add_argument('infile')
    return p


def heatmap(args):

    from seqpy.core.bioio.tabutils import read_file
    from matplotlib import pyplot as plt
    from itertools import combinations

    import seaborn as sns

    df = read_file(args.infile)

    df.set_index('SAMPLE', inplace=True)

    ax = sns.heatmap(df, cmap='gist_heat_r', xticklabels=1, square=True,
                     cbar_kws={"shrink": .8}, yticklabels=1)
    [t.set_fontsize(5) for t in ax.get_yticklabels()]
    [t.set_fontsize(5) for t in ax.get_xticklabels()]
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(args.outplot)



def main(args):
    heatmap(args)


# EOF
