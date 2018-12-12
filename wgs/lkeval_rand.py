# lkeval_rand.py
#
# likelihood evaluator with random SNP selection algorithm


from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from seqpy.core.cfuncs import genoutils, lk_eval
from seqpy.core.bioio import naltparser

import random, itertools

def init_argparser():
    p = arg_parser("Evaluate a set of SNPs for its likelihood-based assignment")
    p.add_argument('-k', default='1000,500,100', help='k list')
    p.add_argument('-o', '--outfile', default='out.accuracy.txt')

    p = naltparser.init_argparser(p)

    return p

def main(args):
    lkeval_rand( args )


def lkeval_rand( args ):

    nalt_parser = naltparser.NAltLineParser( args, datatype='nalt')

    # we need group info
    nalt_parser.parse_grouping()

    group_keys = nalt_parser.group_parser.group_keys

    # remove samples from group which less than certain number
    suitable_groups = set( [ g for g in nalt_parser.group_parser.groups
                                if len(nalt_parser.group_parser.groups[g]) > 2])

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

        result = evaluate_rand( haplotypes, group_keys, k=k, n_inner=10, n_outer=10)
        accuracies.append( lk_eval.calculate_accuracy(result, k) )

    result = pd.concat( accuracies )
    result.to_csv(args.outfile, sep='\t', index=False)


def evaluate_rand( haplotypes, group_keys, k, n_inner=10, n_outer=10 ):
    """ perform iteration to do evaluation """

    group_set = set(group_keys)

    estimation = []

    c = 0
    for i in range(n_outer):
        X_train, X_test, y_train, y_test = train_test_split(haplotypes, group_keys,
                                test_size=0.33, stratify=group_keys)

        #import IPython; IPython.embed()
        for j in range(n_inner):
            c += 1

            # the section below is algorithm-spesific SNP selection

            # we perform randomize here
            L = np.random.randint(0, len(haplotypes[0]), k)
            X_train_0 = X_train[:,L]
            X_test_0 = X_test[:,L]

            # end of algorithm-spesific SNP selection

            evaluator = lk_eval.LikelihoodEvaluator(X_train_0, y_train, X_test_0, y_test)
            evaluator.prepare_matrix()

            estimation.append( evaluator.evaluate() )

            if c % 10 == 0:
                cerr('[I - iteration %d]' % c)


    return estimation







