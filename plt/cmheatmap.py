#!/usr/bin/env spcli
# cmheatmap.py - plot heatmap from confusion matrix

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser("plot heatmap from confusion matrix")
    p.add_argument('--truefile', default=None)
    p.add_argument('--orderfile', default=None,
                   help="File containing lines of class labels.")
    p.add_argument('--annotate', default=False, action='store_true',
                   help="Annotate map with numbers")
    p.add_argument('-o', '--outplot', default='outheatmap.pdf')
    p.add_argument('infile')

    return p


def read_orderfile(orderfile):
    """ return (list of ordered text, set of ordered text) """
    countries = []
    noted_countries = []
    with open(orderfile) as fin:
        for line in fin:
            line = line.strip()
            if line:
                if line.endswith('*'):
                    line = line[:-1]
                    noted_countries.append(line)
                countries.append(line)
    return countries, set(noted_countries)


def reorder_labels(labels, order):
    """ reorder labels based on order list """
    order_set = set(order)
    label_set = set(labels)

    labels_not_in_order_set = label_set - order_set
    ordered_labels = []
    for lbl in order:
        if lbl in label_set:
            ordered_labels.append(lbl)
    ordered_labels += list(labels_not_in_order_set)

    return ordered_labels


def cmheatmap(args):

    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.metrics import confusion_matrix

    # load infile
    true_column = pred_column = None
    infile = args.infile
    if ':' in infile:
        infile, columns = infile.split(':', 1)
        if ',' in columns:
            true_column, pred_column = columns.split(',', 1)
        else:
            pred_column = columns

    truefile = args.truefile
    if truefile is None and pred_column is None:
        cexit('ERR - either set truefile using --truefile or set true column')

    df_pred = pd.read_csv(infile)
    if truefile:
        if ':' in truefile:
            truefile, true_column = truefile.split(':', 1)
        df_true = pd.read_csv(truefile, sep=None)
    else:
        df_true = df_pred

    series_true = df_true[true_column].astype('category')
    series_pred = df_pred[pred_column].astype('category')

    noted_labels = None
    if args.orderfile:
        # reorder label
        order, noted_labels = read_orderfile(args.orderfile)
        ordered_true = reorder_labels(series_true.cat.categories, order)
        series_true = series_true.cat.reorder_categories(ordered_true, ordered=True)
        ordered_pred = reorder_labels(series_pred.cat.categories, order)
        series_pred = series_pred.cat.reorder_categories(ordered_pred, ordered=True)

    sns.set(font_scale=0.6)
    cm = confusion_matrix(series_true.cat.codes, series_pred.cat.codes, normalize='true')

    # modify cm and label order
    n_true = len(series_true.cat.categories)
    n_pred = len(series_pred.cat.categories)

    if n_true < n_pred:
        cm = cm[:n_true, :]
    elif n_pred < n_true:
        cm = cm[:, :n_pred]

    if args.annotate:
        annotate_cm = np.vectorize(lambda v: ('%.3f' % v) if v > 0 else '')
        annotations = annotate_cm(cm)
    else:
        annotations = None

    # plot heatmap
    ax = sns.heatmap(
        cm,
        xticklabels=series_pred.cat.categories,
        yticklabels=series_true.cat.categories,
        annot=annotations,
        annot_kws={'fontsize': 3},
        fmt='',
        cmap='Blues')

    # set ylabel color
    if noted_labels:
        for lbl in ax.get_yticklabels():
            if lbl.get_text() not in noted_labels:
                lbl.set_color('r')

    plt.ylabel('Origin')
    plt.xlabel('Prediction')
    plt.tight_layout()
    plt.savefig(args.outplot)


def main(args):
    cmheatmap(args)

# EOF
