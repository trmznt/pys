""" lkeval_dt.py
    likelihood evaluator using decision tree algorithm
"""
import itertools

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import lk_eval


def init_argparser():
    parser = arg_parser("Evaluate a set of SNPs for its likelihood-based assignment")
    parser.add_argument('-k', default='100,50', help='k list')
    parser.add_argument('-o', '--outfile', default='out.accuracy.dt.txt')
    parser = naltparser.init_argparser(parser)

    return parser


def evaluate_dt(haplotypes, group_keys, k, n_inner=100, n_outer=1):
    """ perform iteration to do evaluation """
    estimation = []
    c = 0
    for _ in range(n_outer):
        X_train, X_test, y_train, y_test = \
                train_test_split(haplotypes, group_keys, test_size=0.33, stratify=group_keys)

        for _ in range(n_inner):
            c += 1

            classifier = DecisionTreeClassifier(class_weight='balanced', max_features=k)
            classifier = classifier.fit(X_train, y_train)
            features = classifier.tree_.feature

            L = np.delete(features, np.where(features == -2))
            X_train_0 = X_train[:, L]
            X_test_0 = X_test[:, L]

            evaluator = lk_eval.LikelihoodEvaluator(X_train_0, y_train, X_test_0, y_test)
            evaluator.prepare_matrix()

            estimation.append(evaluator.evaluate())

            if c % 10 == 0:
                cerr('[I - iteration %d]' % c)

    return estimation


def lkeval_dt(args):
    nalt_parser = naltparser.NAltLineParser(args, datatype='nalt')

    # we need group info
    nalt_parser.parse_grouping()
    group_keys = nalt_parser.group_parser.group_keys

    # remove samples from group which less than certain number
    suitable_groups = {g for g in nalt_parser.group_parser.groups
                       if len(nalt_parser.group_parser.groups[g]) > 2}

    cerr('[I - reading %s samples]' % len(nalt_parser.samples))
    mask = [False] * len(nalt_parser.samples)
    for i in range(len(nalt_parser.samples)):
        if group_keys[i] in suitable_groups:
            mask[i] = True

    samples = list(itertools.compress(nalt_parser.samples, mask))
    group_keys = list(itertools.compress(group_keys, mask))
    cerr('[I - masking to %d samples]' % len(samples))

    region = nalt_parser.parse_whole(mask=mask)
    haplotypes = region.haplotypes()

    accuracies = []
    for k in [ int(x) for x in args.k.split(',') ]:

        cerr('[I - iterating for k = %d]' % k)
        result = evaluate_dt( haplotypes, group_keys, k=k, n_inner=10, n_outer=10)
        accuracies.append( lk_eval.calculate_accuracy(result, k) )

    result = pd.concat( accuracies )
    result.to_csv(args.outfile, sep='\t', index=False)


def main(args):
    lkeval_dt(args)
