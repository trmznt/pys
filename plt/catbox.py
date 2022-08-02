#!/usr/bin/env spcli
"""
    plot_accuracies.py

"""

import sys
import pandas as pd
import seaborn as sns

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    parser = arg_parser("plot_accuracies.py")
    parser.add_argument('--variables', default='MCC')
    parser.add_argument('--labelcolumn', default='CLASS')
    parser.add_argument('--labels', default=None)
    parser.add_argument('--hue', default=None)
    parser.add_argument('--method', default='')
    parser.add_argument('--modelid', default='')
    parser.add_argument('--region', default=None)
    parser.add_argument('--row', default=None)
    parser.add_argument('--aspect', type=int, default=3)
    parser.add_argument('--height', type=float, default=-1)
    parser.add_argument('--order', default=None, help="File containing order of x axis")
    parser.add_argument('--statsummary', default='', help='Add stat summary, group by this column')
    parser.add_argument('--outplot', default='outplot.pdf')
    parser.add_argument('--palette', default='muted',
                        help='Available colour palette: Paired pastel')
    parser.add_argument('infile', nargs='+')

    return parser


def drop_dataframe(dataframe, variables, label_column, labels=None):
    # untouchable = ['REG', '_k']
    # untouchable.extend(variables)
    # for column in dataframe.columns:
    #    if column in untouchable:
    #        continue
    #    dataframe.drop(columns=column, inplace=True)

    df = dataframe[variables + [label_column]]

    if labels is not None:
        df = df[df[label_column].isin(labels)]

    return df


def calculate_required_aspect(data, label_column):
    total_width = 1.37 * len(data[label_column].unique())
    total_height = total_width / 1.42 / 2   # bad approximation of square root of 2
    cerr(f'[I - Using figure size: {total_height:.2f} x {total_width:.2f} inches]')

    return total_width, total_height


def read_infiles(args):

    dataframes = []
    for infile in args.infile:

        df = pd.read_csv(infile, header=0, sep='\t')

        # modify header
        if 'k' in df.columns:
            df.rename(columns={"k": "label"}, inplace=True)
        if '_k' in df.columns:
            df.rename(columns={'_k': 'k'}, inplace=True)

        dataframes.append(df)

    if len(dataframes) == 1:
        return dataframes[0]

    return pd.concat(dataframes)


def plot_accuracies(args):

    data = read_infiles(args)

    if args.statsummary:

        cerr(f'[Adding stats grouped by {args.statsummary}]')

        group_index = [args.statsummary, args.hue] if args.hue else [args.statsummary]
        group_columns = group_index + [args.variables]

        df = data[group_columns].groupby(group_index).agg(['median', 'min']).stack().reset_index()
        df.rename(columns={'level_2': args.labelcolumn}, inplace=True)

        data = pd.concat([data, df], ignore_index=True)

    # filter data
    if args.method:
        data = data[data['METHOD'] == args.method]

    if args.modelid:
        modelids = args.modelid.split(',')
        data = data[data['MODELID'].isin(modelids)]

    # parse variables & regions

    variables = args.variables.split(',') + ([args.hue] if args.hue else [])
    if args.row:
        variables += [args.row]

    if args.labels:
        labels = args.labels.split(',')
    else:
        labels = None

    if args.order:
        order = open(args.order).read().split('\n')
        if args.statsummary:
            order += ['median', 'min']
    else:
        order = None

    if args.hue:
        melted = [args.labelcolumn, args.hue]
        if args.row:
            melted += [args.row]
        cerr(f'[I - melted variables: {melted}')
        data = drop_dataframe(data, variables, args.labelcolumn, labels).melt(melted)

    sns.set_palette(args.palette)

    if labels is None:
        if args.height > 0:
            height = args.height
        else:
            width, height = calculate_required_aspect(data, args.labelcolumn)
        row = args.row if args.row else 'variable'
        plots = sns.catplot(x=args.labelcolumn, y='value', data=data, hue=args.hue, row=row,
                            kind='box', sharex=True, height=height, aspect=args.aspect,
                            order=order, showfliers=False, legend_out=True, linewidth=1)
    else:
        plots = sns.catplot(x=args.labelcolumn, y='value', data=data, hue=args.hue, row='variable',
                            kind='box', sharex=True)
    axes = plots.axes
    for axis in axes:
        axis = axis[0]
        title = (axis.get_title()).split(' = ')[-1]
        axis.set_title('')
        axis.set_ylabel(title, fontsize='x-large')

    # rotate x ticks
    [(x.set_rotation(60), x.set_fontsize(x.get_fontsize() * 1.5)) for x in plots.axes[-1][0].get_xticklabels()]
    plots.axes[-1][0].set_xlabel('')

    fig = plots.fig

    # import IPython; IPython.embed()

    # this will ensure that labels for x-axis are shown, but currently reposition the legend as well
    # fig.tight_layout(w_pad=5)
    fig.savefig(args.outplot, dpi=300, bbox_inches='tight')


def main(args):

    plot_accuracies(args)

# EOF
