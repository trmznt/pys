#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


# note on matplotlib default category10 palette (from Vega and d3)
# ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
# '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outplot')
    p.add_argument('-c', '--colorfile',
                   help=("a TSV/CSV file containg the group, COLOR and MARKER columns "
                         "for defining the color and marker of each data point"))
    p.add_argument('--hue', default=None,
                   help="Column name to used to differentiate data points")
    p.add_argument('--dpi', type=int, default=600)
    p.add_argument('--jitter', type=float, default=0.005)
    p.add_argument('--fontscale', type=float, default=1,
                   help="font scale, default=1")
    p.add_argument('-s', '--size', type=int, default=10,
                   help="size of markers")
    p.add_argument('-w', '--width', type=float, default=8.0,
                   help="Width & height of each scatter plot, in inch. Default is 8.")
    p.add_argument('--columns', required=True,
                   help="Column names to be used as x-axis and y-axis data")
    p.add_argument('--legend_out', action='store_true', default=False,
                   help="Legend will be put outside (right of) plot")
    p.add_argument('infile')
    return p


def scatter(args):

    from seqpy.core.bioio.tabutils import read_file
    from matplotlib import pyplot as plt
    from itertools import combinations

    import seaborn as sns

    df = read_file(args.infile)
    columns = args.columns.split(',')

    # setting up color and markers

    colors = None
    markers = True
    hue_order = None

    if args.colorfile:
        if args.hue is None:
            cexit('[Using colorfile requires --hue option as well]')
        color_df = read_file(args.colorfile)
        if 'COLOR' in color_df.columns:
            colors = color_df.meta.to_dict(args.hue, 'COLOR')
        if 'MARKER' in color_df.columns:
            cerr('Markers will be set based on column MARKER')
            markers = color_df.meta.to_dict(args.hue, 'MARKER')
        hue_order = color_df[args.hue]

    data_axes = list(combinations(columns, 2))
    N = len(data_axes)

    sns.set_theme(font_scale=args.fontscale, style='ticks')

    fig = plt.figure(figsize=(args.width + (2.5 if args.legend_out else 0),
                              args.width * N), dpi=args.dpi)
    fig_idx = 1

    for x_ax, y_ax in data_axes:

        x_jitter = max(df[x_ax]) - min(df[x_ax]) * args.jitter
        y_jitter = max(df[y_ax]) - min(df[y_ax]) * args.jitter

        ax = fig.add_subplot(N, 1, fig_idx)
        fig_idx += 1

        sns.scatterplot(x=x_ax, y=y_ax, data=df,
                        hue=args.hue, style=args.hue, hue_order=hue_order,
                        palette=colors,
                        markers=markers,
                        s=args.size,
                        x_jitter=x_jitter,
                        y_jitter=y_jitter,
                        alpha=0.9,
                        )
        if args.legend_out:
            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    fig.tight_layout()
    fig.savefig(args.outplot)
    cerr(f'[Scatter plot written to {args.outplot}]')


def main(args):
    scatter(args)

# EOF
