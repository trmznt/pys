#!/usr/bin/env spcli
"""
    plot_accuracies.py

"""

import sys

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
    parser.add_argument('--palette', default='muted',
                        help='Available colour palette: Paired pastel')
    parser.add_argument('--rotation', type=int, default=60,
                        help='Label rotation in degree (default: 60)')

    parser.add_argument('-i', '--interactive', default=False, action='store_true',
                        help='Drop to IPython interactive session after before exit.')

    parser.add_argument('--outtable', default='')
    parser.add_argument('--outplot', default='outplot.pdf')

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

    import pandas as pd

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

    import pandas as pd
    import seaborn as sns

    data = read_infiles(args)

    if args.outtable:

        # write numbers to tables

        stacked_column = 'level_2'
        group_index = [args.labelcolumn, args.hue] if args.hue else [args.labelcolumn]
        group_columns = group_index + [args.variables]

        if args.row:
            group_index.append(args.row)
            group_columns.append(args.row)
            stacked_column = 'level_3'

        df = data[group_columns].groupby(group_index).agg(['median', 'min']).stack().reset_index()
        # df now is:
        #            CLASS       MODEL level_2       MCC
        # 0    Afghanistan        BR38  median  0.834963
        # 1    Afghanistan        BR38     min  0.361357
        # 2    Afghanistan       GEO33  median  0.849747
        # 3    Afghanistan       GEO33     min  0.579221
        # 4    Afghanistan  GEO33+BR38  median  0.946015
        # ..           ...         ...     ...       ...
        # 289      Vietnam  GEO50+BR38     min  0.212224
        # 290      Vietnam       GEO55  median  0.756743
        # 291      Vietnam       GEO55     min  0.277208
        # 292      Vietnam  GEO55+BR38  median  0.756743
        # 293      Vietnam  GEO55+BR38     min  0.429472

        # df.rename(columns={'level_2': args.labelcolumn}, inplace=True)

        if args.row:
            for val in df[args.row].unique():
                df2 = df[df[args.row] == val]
                df2_med = df2[df2[stacked_column] == 'median'].pivot(index='CLASS',
                                                                     columns='MODEL',
                                                                     values='MCC')
                df2_min = df2[df2[stacked_column] == 'min'].pivot(index='CLASS',
                                                                  columns='MODEL',
                                                                  values='MCC')

                df2_med.reset_index().to_csv(f'{args.outtable}-{val}-med.tsv',
                                             index=False, sep='\t')
                df2_min.reset_index().to_csv(f'{args.outtable}-{val}-min.tsv',
                                             index=False, sep='\t')

        else:
            df_med = df[df[stacked_column] == 'median'].pivot(index='CLASS',
                                                              columns='MODEL',
                                                              values='MCC')
            df_min = df[df['level_2'] == 'min'].pivot(index='CLASS',
                                                      columns='MODEL',
                                                      values='MCC')

            df_med.reset_index().to_csv(args.outtable + '-med.tsv', index=False, sep='\t')
            df_min.reset_index().to_csv(args.outtable + '-min.tsv', index=False, sep='\t')

    if args.statsummary:

        cerr(f'[Adding stats grouped by {args.statsummary}]')

        stacked_column = 'level_2'
        group_index = [args.statsummary, args.hue] if args.hue else [args.statsummary]
        group_columns = group_index + [args.variables]

        if args.row:
            group_index.append(args.row)
            group_columns.append(args.row)
            stacked_column = 'level_3'

        df = data[group_columns].groupby(group_index).agg(['median', 'min']).stack().reset_index()
        df.rename(columns={stacked_column: args.labelcolumn}, inplace=True)

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
    [(x.set_rotation(args.rotation), x.set_fontsize(x.get_fontsize() * 1.5))
        for x in plots.axes[-1][0].get_xticklabels()]
    plots.axes[-1][0].set_xticklabels(
        plots.axes[-1][0].get_xticklabels(),
        horizontalalignment='right')
    plots.axes[-1][0].set_xlabel('')

    fig = plots.fig

    # import IPython; IPython.embed()

    # this will ensure that labels for x-axis are shown, but currently reposition the legend as well
    # fig.tight_layout(w_pad=5)
    fig.savefig(args.outplot, dpi=300, bbox_inches='tight')

    if args.interactive:
        import IPython
        IPython.embed()


def main(args):

    plot_accuracies(args)

# EOF
