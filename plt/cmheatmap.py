#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

def init_argparser():
    p = arg_parser("plot heatmap from confusion matrix")
    p.add_argument('--truefile', default=None)
    p.add_argument('-o', '--outplot', default='outheatmap.pdf')
    p.add_argument('infile')

    return p

def cmheatmap( args ):

    import pandas as pd, seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.metrics import confusion_matrix

    # load infile
    true_column = pred_column = None
    infile = args.infile
    if ':' in infile:
        infile, columns = infile.split(':', 1)
        if ',' in columns:
            true_column, pred_column = columns.split(',',1)
        else:
            pred_column = columns

    truefile = args.truefile
    if truefile is None and pred_column is None:
        cexit('ERR - either set truefile using --truefile or set true column')

    if ':' in truefile:
        truefile, true_column = truefile.split(':', 1)

    df_pred = pd.read_csv(infile, sep='\t')
    if truefile:
        df_true = pd.read_csv(truefile, sep='\t')
    else:
        df_true = df_pred

    # recoding
    #series_true = df_true[true_column].astype('category')
    #series_pred = df_pred[pred_column].astype('category')
    #cm = confusion_matrix(series_true.cat.codes, series_pred.cat.codes, normalize='true')
    cm = confusion_matrix(df_true[true_column], df_pred[pred_column], normalize='true')
    #import IPython; IPython.embed()
    sns.heatmap(cm,
                xticklabels = series_pred.cat.categories,
                yticklabels = series_true.cat.categories,
                cmap='Blues')

    plt.ylabel('Origin')
    plt.xlabel('Prediction')
    plt.tight_layout()
    plt.savefig(args.outplot)

def main( args ):
    cmheatmap( args )