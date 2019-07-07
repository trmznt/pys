"""
    plot_accuracies.py

"""

import sys
import pandas as pd, seaborn as sns

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

def init_argparser():
    parser = arg_parser("plot_accuracies.py")
    parser.add_argument('--variables', default='F')
    parser.add_argument('--method', default='')
    parser.add_argument('--region', default=None)
    parser.add_argument('--hue', default='k')
    parser.add_argument('--aspect', type=int, default=3)
    parser.add_argument('--order', default=None, help="File containing order of x axis")
    parser.add_argument('--outplot', default='outplot.pdf')
    parser.add_argument('infile', nargs='+')

    return parser


def drop_dataframe(dataframe, variables, regions=None):
    #untouchable = ['REG', '_k']
    #untouchable.extend(variables)
    #for column in dataframe.columns:
    #    if column in untouchable:
    #        continue
    #    dataframe.drop(columns=column, inplace=True)

    df = dataframe[ variables + ['REG', 'k'] ]

    if regions is not None:
        df = df[df['REG'].isin(regions)]

    return df


def calculate_required_aspect(data):
    total_width = 1.37 * len(data['REG'].unique())
    total_height = total_width / 1.42 / 2 # bad approximation of square root of 2
    print('Using figure size: {:.2f} x {:.2f} inches'.format(
              total_height, total_width),
          file=sys.stderr)

    return total_width, total_height


def read_infiles(args):

    dataframes = []
    for infile in args.infile:

        df = pd.read_csv(infile, header=0, sep='\t')

        # modify header
        if 'k' in df.columns:
            df.rename( columns= {"k": "label"}, inplace=True)
        if '_k' in df.columns:
            df.rename( columns= {'_k': 'k'}, inplace=True)

        dataframes.append( df )

    if len(dataframes) == 1:
        return dataframes[0]

    return pd.concat( dataframes )

def plot_accuracies(args):


    data = read_infiles(args)


    # filter data
    if args.method:
        data = data[ data['METHOD'] == args.method ]

    # parse variables & regions

    variables = args.variables.split(',') + [ args.hue ]
    print(variables)

    if args.region:
        regions = args.region.split(',')
    else:
        regions = None

    if args.order:
        order = open(args.order).read().split('\n')
    else:
        order = None

    data = drop_dataframe(data, variables, regions).melt(['REG', 'k', args.hue])

    sns.set_palette('Paired')

    if regions is None:
        width, height = calculate_required_aspect(data)
        plots = sns.catplot(x='REG', y='value', data=data, hue=args.hue, row='variable',
                            kind='box', sharex=True, height=height, aspect=args.aspect,
                            order=order, showfliers=False, legend_out=True)
    else:
        plots = sns.catplot(x='REG', y='value', data=data, hue=args.hue, row='variable',
                            kind='box', sharex=True)
    axes = plots.axes
    for axis in axes:
        axis = axis[0]
        title = (axis.get_title()).split(' = ')[-1]
        axis.set_title('')
        axis.set_ylabel(title, fontsize='x-large')

    # rotate x ticks
    [ (x.set_rotation(60), x.set_fontsize( x.get_fontsize() * 2)) for x in plots.axes[-1][0].get_xticklabels() ]
    plots.axes[-1][0].set_xlabel('')

    fig = plots.fig

    #import IPython; IPython.embed()

    # this will ensure that labels for x-axis are shown, but currently reposition the legend as well
    #fig.tight_layout(w_pad=5)
    fig.savefig( args.outplot, dpi=300, bbox_inches='tight')


def main(args):

    plot_accuracies(args)
