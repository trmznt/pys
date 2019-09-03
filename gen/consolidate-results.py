#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser
from seqpy.core.cfuncs import lkmodels

import numpy as np
import pandas as pd
import pickle

def init_argparser(p=None):

    if p is None:
        p = arg_parser('consolidate_predictions')

    p = grpparser.init_argparser(p)
    p.add_argument('--fmt', default='text', choices=['pickle', 'npy', 'text', 'list'])
    p.add_argument('--samplefile', default='')
    p.add_argument('--threshold', type=float, default=1.0)
    p.add_argument('--outreport', default='')
    p.add_argument('infile')
    return p


def main( args ):
    consolidate_predictions( args )


def consolidate_predictions( args ):

    outreport = None

    if args.samplefile:
        samples = read_samplefile( args.samplefile, args.fmt)
    else:
        samples = None

    group_parser = grpparser.GroupParser( args )
    group_parser.assign_groups(samples)
    #group_parser.group_keys contains [ 'grp1', 'grp2', etc]
    group_keys = group_parser.group_keys

    with open(args.infile, 'rb') as f:
        predictions = pickle.load(f)

    if args.outreport:
        outreport = open(args.outreport, 'wb')
        from sklearn.metrics import confusion_matrix
    reports = {}

    normalize = True

    for model in predictions:
        model_pred = predictions[model]

        for k in model_pred:

            cerr('Preparing for model: {} k: {}'.format(model, k))
            df = generate_dataframe( model_pred[k])

            group_indexes = np.argmax(df.values, axis=1)
            group_predictions = df.columns[ group_indexes[:,None] ]
            for i in range(len(group_indexes)):
                predicted_group = df.columns[group_indexes[i]]
                prediction_confidence = df.values[i, group_indexes[i]]
                if prediction_confidence < args.threshold or predicted_group != group_keys[i]:
                    cout('{}: {} -> {} ({})'.format(samples[i], group_keys[i], predicted_group, prediction_confidence))

            if outreport:

                score = lkmodels.calculate_scores(group_keys, group_predictions)
                confmat = confusion_matrix(group_keys, group_predictions)

                if normalize:
                    confmat = confmat.astype('float') / confmat.sum(axis=1)[:, np.newaxis]
                    cerr("[I - Normalized confusion matrix]")
                else:
                    cerr('[I - Confusion matrix, without normalization]')

                reports['{}|{}'.format(model, k)] = {
                    'score': score
                    , 'confmat': confmat
                }

    if outreport:
        pickle.dump(reports, outreport)
        cerr('[I - writing pickled report to {}]'.format(args.outreport))




def read_samplefile( infile, fmt ):

    if fmt in ['pickle', 'npy']:
        from seqpy.core.bioio import naltparser
        from types import SimpleNamespace
        nalt_args = SimpleNamespace( infile = infile, fmt = fmt, n = -1)
        nalt_parser = naltparser.NAltLineParser(nalt_args
                , with_group=False, with_position=False)
        samples = nalt_parser.samples

    elif fmt == 'list':
        with gzopen(infile) as f:
            buf = f.read()
            samples = buf.strip().split()

    else:
        # only read the first line
        with gzopen(infile) as f:
            samples = f.readline().strip().split()

    return samples


def generate_dataframe( pred_list ):

    labels = np.unique( pred_list )
    df = pd.DataFrame( np.zeros( (len(pred_list[0]), len(labels)) ), columns = labels)

    for a_list in pred_list:
        for i in range(len(a_list)): df.loc[i, a_list[i]] += 1

    # convert from count to frequency
    totals = df.values.sum(axis=1)
    values = df.values/totals[:,None]

    return pd.DataFrame( values, columns = labels )
