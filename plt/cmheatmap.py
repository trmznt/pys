#!/usr/bin/env spcli
# cmheatmap.py - plot heatmap from confusion matrix

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np
import pandas as pd


def init_argparser():
    p = arg_parser("plot heatmap from confusion matrix")
    p.add_argument('--truefile', default=None)
    p.add_argument('--orderfile', default=None,
                   help="File containing lines of class labels.")
    p.add_argument('--annotate', default=False, action='store_true',
                   help="Annotate map with numbers")
    p.add_argument('--colorlabelaxis', default='')
    p.add_argument('--labelfontscale', type=float, default=1.25)
    p.add_argument('--fontscale', type=float, default=1)
    p.add_argument('--annotfontscale', type=float, default=0.5)
    p.add_argument('--titlefontscale', type=float, default=1.75)
    p.add_argument('--titlexpos', type=float, default=-0.2)
    p.add_argument('--nocompact', default=False, action='store_true')
    p.add_argument('--outtable', default='')
    p.add_argument('--title', default='')
    p.add_argument('-o', '--outplot', default='outheatmap.pdf')
    p.add_argument('infile',
                   help="prediction result file, can be supplemented with column names, "
                   "eg: filename:PREDCOLUMN or filename:TRUECOLUMN,PREDCOLUMN")

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


def prepare_dataframe_metrics(y_true, y_pred, labels):

    from sklearn.metrics._classification import multilabel_confusion_matrix
    from sklearn.utils.multiclass import unique_labels

    labels = unique_labels(y_true)

    MCM = multilabel_confusion_matrix(
        y_true,
        y_pred,
        labels=labels
    )

    tn_sum = MCM[:, 0, 0]
    fn_sum = MCM[:, 1, 0]
    tp_sum = MCM[:, 1, 1]
    fp_sum = MCM[:, 0, 1]
    pp_sum = tp_sum + fp_sum
    pn_sum = fn_sum + tn_sum
    p_sum = tp_sum + fn_sum
    n_sum = fp_sum + tn_sum

    # MCC
    denom = np.sqrt((tp_sum + fp_sum) * (tp_sum + fn_sum) * (tn_sum + fp_sum) * (tn_sum + fn_sum))
    denom[denom == 0.0] = 1  # avoid division by 0
    mcc = (tp_sum * tn_sum - fp_sum * fn_sum) / denom

    return pd.DataFrame({'CLASS': labels,
                         'TN': tn_sum, 'FN': fn_sum, 'TP': tp_sum, 'FP': fp_sum,
                         'PP': pp_sum, 'PN': pn_sum, 'P': p_sum, 'N': n_sum,
                         'MCC': mcc}), MCM


def cmheatmap(args):

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

    df_pred = pd.read_csv(infile, sep=None)
    if truefile:
        if ':' in truefile:
            truefile, true_column = truefile.split(':', 1)
        df_true = pd.read_csv(truefile, sep=None)
    else:
        df_true = df_pred

    series_true = df_true[true_column].astype('category')
    series_pred = df_pred[pred_column].astype('category')

    noted_labels = order = None
    if args.orderfile:
        # reorder label
        order, noted_labels = read_orderfile(args.orderfile)
        ordered_true = ordered_pred = order
        if not args.nocompact:
            ordered_true = reorder_labels(series_true.cat.categories, order)
            ordered_pred = reorder_labels(series_pred.cat.categories, order)
        #import IPython; IPython.embed()
        series_true = series_true.cat.set_categories(ordered_true, ordered=True)
        series_pred = series_pred.cat.set_categories(ordered_pred, ordered=True)

    if args.nocompact:

        cm = confusion_matrix(series_true, series_pred, labels=order, normalize='true')

    else:

        cm = confusion_matrix(series_true.cat.codes, series_pred.cat.codes, normalize='true')

        # modify cm and label order
        n_true = len(series_true.cat.categories)
        n_pred = len(series_pred.cat.categories)

        if n_true < n_pred:
            cm = cm[:n_true, :]
        elif n_pred < n_true:
            cm = cm[:, :n_pred]

    if args.annotate:
        annotate_cm = np.vectorize(lambda v: ('%.3f' % v) if 0 < v < 0.5 else '')
        annotations = annotate_cm(cm)
    else:
        annotations = None

    # plot heatmap
    sns.set(font_scale=args.fontscale)
    fig, ax = plt.subplots(1)
    fontsize = ax.xaxis.label.get_fontsize()

    ax = sns.heatmap(
        cm,
        xticklabels=series_pred.cat.categories,
        yticklabels=series_true.cat.categories,
        annot=annotations,
        annot_kws={'size': fontsize * args.annotfontscale},
        fmt='',
        cmap='Blues',
        ax=ax
    )

    # set ylabel color
    if noted_labels:
        for ca in args.colorlabelaxis:
            func = ax.get_yticklabels if ca == 'y' else ax.get_xticklabels
            for lbl in func():
                if lbl.get_text() not in noted_labels:
                    lbl.set_color('r')

    [(x.set_rotation(90), x.set_fontsize(fontsize))
        for x in ax.get_xticklabels()]

    [(y.set_rotation(0), y.set_fontsize(fontsize))
        for y in ax.get_yticklabels()]

    ax.set_ylabel('Origin', fontsize=fontsize * args.labelfontscale)
    ax.set_xlabel('Prediction', fontsize=fontsize * args.labelfontscale)
    ax.set_title(args.title, fontsize=fontsize * args.titlefontscale, x=args.titlexpos)
    fig.tight_layout()
    fig.savefig(args.outplot)

    if args.outtable:

        df, mcm = prepare_dataframe_metrics(series_true, series_pred, None)
        df.to_csv(args.outtable, index=False, sep='\t')


def main(args):
    cmheatmap(args)

# EOF
