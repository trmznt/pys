#!/usr/bin/env spcli
# cmheatmap.py - plot heatmap from confusion matrix

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np
import pandas as pd

from itertools import zip_longest

import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

import math


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
    p.add_argument('--colorbarscale', type=float, default=1.0)
    p.add_argument('--titlexpos', type=float, default=-0.2)
    p.add_argument('--compactaxis', default='')
    p.add_argument('--xrotation', type=float, default=45)
    p.add_argument('--width', type=float, default=15)
    p.add_argument('--height', type=float, default=20)
    p.add_argument('--outtable', default='')
    p.add_argument('--title', default='')
    p.add_argument('-o', '--outplot', default='outheatmap.pdf')
    p.add_argument('infiles', nargs='+',
                   help="prediction result file, can be supplemented with column names, "
                   "eg: filename:PREDCOLUMN or filename:TRUECOLUMN,PREDCOLUMN")

    return p


def read_orderfile(orderfile):
    """ return (list of ordered text, set of ordered text) """
    labels = []
    noted_labels = []
    with open(orderfile) as fin:
        for line in fin:
            line = line.strip()
            if line:
                if line.endswith('*'):
                    line = line[:-1]
                    noted_labels.append(line)
                labels.append(line)
    return np.array(labels), set(noted_labels)


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


def cmheatmap(args, ax, infile, title, index, axbar=None):

    # load infile
    true_column = pred_column = None
    if ':' in infile:
        infile, columns = infile.split(':', 1)
        if ',' in columns:
            true_column, pred_column = columns.split(',', 1)
        else:
            pred_column = columns

    truefile = args.truefile
    if truefile is None and pred_column is None:
        cexit('ERR - either set truefile using --truefile or set true column')

    df_pred = pd.read_csv(infile, sep=None, engine='python')
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
    else:
        order = np.array(sorted(set(series_true) + set(series_pred)))

    cm = confusion_matrix(series_true, series_pred, labels=order, normalize='true')

    xlabels = ylabels = order

    if 'x' in args.compactaxis:
        idx_x = cm.sum(axis=0) > 0
        cm = cm[:, idx_x]
        xlabels = order[idx_x]

    if 'y' in args.compactaxis:
        #import IPython; IPython.embed()
        idx_y = cm.sum(axis=1) > 0
        cm = cm[idx_y, :]
        ylabels = order[idx_y]

    if args.annotate:
        annotate_cm = np.vectorize(lambda v: ('%.3f' % v) if 0 < v < 2 else '')
        annotations = annotate_cm(cm)
    else:
        annotations = None

    # plot heatmap
    fontsize = ax.xaxis.label.get_fontsize()

    cbar = False
    cbar_ax = None
    if axbar is not None:
        cerr('Adding color bar')
        cbar = True
        cbar_ax = axbar

    ax = sns.heatmap(
        cm,
        xticklabels=xlabels,
        yticklabels=ylabels,
        annot=annotations,
        annot_kws={'size': fontsize * args.annotfontscale},
        fmt='',
        cmap='Blues',
        square=True,
        cbar_kws=dict(shrink=args.colorbarscale),
        ax=ax,
        #cbar=False if (index + 1) % 2 == 0 else True,
        cbar=cbar,
        cbar_ax=cbar_ax,
        linewidths=2,
        linecolor='white',
    )

    ax.tick_params(labelbottom=True, labelleft=True,
                   top=False, bottom=True, left=True, right=False)

    # set ylabel color
    if noted_labels:
        for ca in args.colorlabelaxis:
            func = ax.get_yticklabels if ca == 'y' else ax.get_xticklabels
            for lbl in func():
                if lbl.get_text() not in noted_labels:
                    lbl.set_color('r')

    ax.set_xticklabels(ax.get_xticklabels(), rotation=args.xrotation,
                       fontsize=fontsize, horizontalalignment='right')

    [(y.set_rotation(0), y.set_fontsize(fontsize))
        for y in ax.get_yticklabels()]

    ax.set_xlabel('Prediction', fontsize=fontsize * args.labelfontscale)
    ax.set_ylabel('Origin', fontsize=fontsize * args.labelfontscale)
    ax.set_title(title, fontsize=fontsize * args.titlefontscale, x=args.titlexpos, loc='left')

    if args.outtable:

        df, mcm = prepare_dataframe_metrics(series_true, series_pred, None)
        df.to_csv(f'{args.outtable}-{index:02d}.tsv', index=False, sep='\t')


def multicmheatmap(args):

    # no of infiles
    input_len = len(args.infiles)
    if input_len > 1:
        rows = math.ceil(input_len / 2)
        cols = 2
    else:
        rows = cols = 1

    sns.set(font_scale=args.fontscale)
    fig, axs = plt.subplots(rows, cols, squeeze=False, figsize=(args.width, args.height),
                            sharex=False, sharey=False,)
                            #gridspec_kw=dict(hspace=-0.1))
    axs = axs.flatten()
    #cbar_ax = fig.add_axes([.91, .3, .03, .4])

    for (idx, (ax, infile, title)) in enumerate(zip_longest(axs, args.infiles, args.title.split(';'))):
        labeltop = labelbottom = labelright = labelleft = False
        if not infile:
            #ax.remove()
            continue
        labelleft = labelbottom = True

        cmheatmap(args, ax, infile, title, idx, axbar=axs[-1] if idx == 0  else None)

    axs[-1].set_aspect(10)
    fig.tight_layout()
    fig.savefig(args.outplot)


def main(args):
    multicmheatmap(args)

# EOF
