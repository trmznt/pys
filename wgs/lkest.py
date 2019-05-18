"""
    lkest.py
    likelihood estimator
"""

import itertools, time, math, os, pickle, datetime
from multiprocessing import Pool, RawArray

import numpy as np, pandas as pd, pyparsing as pp
from sklearn.model_selection import StratifiedKFold
from sklearn.tree import DecisionTreeClassifier
import allel

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import lkprof

def init_argparser(p = None):
    if p is None:
        p = arg_parser('lkest.py - SNP-based likelihood estimator')

    p = naltparser.init_argparser(p)

    p.add_argument('-c', '--config', default=None, help = 'config file to use (optional)')
    p.add_argument('--mode', default='predict', choices = mode.keys()
                , help = 'run mode')
    p.add_argument('--logfile', default=None)
    p.add_argument('-o', '--outfile', default='outfile.txt')

    # cross validation params
    p.add_argument('--repeats', default=10, type=int, help='number of repeats')
    p.add_argument('--fold', default=5, type=int, help='number of fold')
    p.add_argument('--iter', default=1, type=int, help='iteration per fold')

    # simulate
    p.add_argument('-j', default=1, type=int, help='no of process')
    p.add_argument('--outsnp', default=None, help='output filename for index of SNPs')
    p.add_argument('-m', '--method', choices=['fix', 'rand', 'list', 'dt', 'hfst', 'hfst+dt'])
    p.add_argument('-k', default='100,50'
            , help='list of k, with k interpreted depending on selected method')
    p.add_argument('--insnp', default=None, help='input file containing index of SNPs')
    p.add_argument('--guidetree', default=None, help='guide tree used by hfst and hfst+dt')
    p.add_argument('--minfst', default=None, help='minimum FST for hfst and hfst+dt')
    p.add_argument('--mac', default=0, help='filtering for MAC')


    # train
    p.add_argument('--outprofile')

    # predict
    p.add_argument('--inprofile')

    return p



class RandomSelector(object):

    code = 'rand'

    def __init__(self, seed=None):
        self.randomstate = np.random.RandomState(seed)
        self.logs = []    # put logs here

    def reseed(self, seed):
        self.randomstate = np.random.RandomState(seed)

    def select(self, haplotypes, groups, haplotest, k):
        """ return (index list of selected SNPs, prediction result) """

        return (self.randomstate.randint(0, len(haplotypes[0]), k), None)

    def log(self, logmsg):
        self.logs.append(logmsg)

    def get_loglines(self):
        return self.logs


class FixSNPSelector(RandomSelector):

    def __init__(self, L=None, snpfile=None):
        if L is not None:
            self.L = L
        elif snpfile:
            self.L = np.array( [ int(x) for x in open(snpfile).read().split('\n')] )
        super().__init__()

    def select(self, haplotypes, groups, haplotest, k=None):
        return (self.L, None)


class DecisionTreeSelector(RandomSelector):

    code = 'dt'

    def select(self, haplotypes, groups, haplotest, k=None):

        if k <= 0:
            k = None

        classifier = DecisionTreeClassifier(class_weight='balanced', max_features=k,
                        random_state = self.randomstate)
        classifier = classifier.fit(haplotypes, groups)
        features = classifier.tree_.feature

        return (np.unique(np.delete(features, np.where(features == -2))),
                    classifier.predict(haplotest))


def parse_guide_tree( treefile ):

    pp.ParserElement.inlineLiteralsUsing(pp.Suppress)
    identity = pp.Word(pp.alphas + ' ')
    element = pp.Group(identity)
    parser = pp.OneOrMore( pp.nestedExpr(content=pp.delimitedList(element) + pp.Optional(','))
                      | pp.delimitedList(element))

    return parser.parseString( treefile.read() ).asList()


def flatten(lst):
    return sum( ([x] if not isinstance(x, list) else flatten(x)
             for x in lst), [] )


def traverse(tree, level=0):

    #print(tree)
    if len(tree) != 2:
        raise RuntimeError('[E - FATAL PROG ERR: misformat tree structure: child size=%d]' % len(tree))

    if isinstance(tree[0], list) and len(tree[0]) > 1:
        pop1 = flatten(tree[0])
        next1 = traverse(tree[0], level+1)
    else:
        pop1 = tree[0]
        next1 = []

    if isinstance(tree[1], list) and len(tree[1]) > 1:
        pop2 = flatten(tree[1])
        next2 = traverse(tree[1], level+1)
    else:
        pop2 = tree[1]
        next2 = []

    return [ (level, pop1, pop2) ] + next1 + next2

def count_allele(haplotypes):
    """ return numpy array of [ [c1, c2], [c1, c2], [c1, c2], .... ] based on haplotypes """

    n_alt = haplotypes.sum(axis=0)
    ac = np.empty( (len(n_alt), 2))
    ac[:,0] = len(haplotypes)*2 - n_alt
    ac[:,1] = n_alt
    return ac


class HierarchicalFSTSelector(RandomSelector):

    code = 'hfst'

    def __init__(self, seed=None, guide_tree=None, min_fst = 0.9):
        super().__init__(seed)
        if guide_tree is None:
            cexit('[E - HierarchicalFSTSelector requires guide tree]')
        self.guide_tree = guide_tree
        self.min_fst = min_fst


    def select(self, haplotypes, groups, haplotest, k=None):

        # we use k for redundancy parameters
        if k == 0 or k is None:
            k = 1

        candidate_L = []     # [ (pos, rank, no_actual_pops)]
        # we traverse through the tree
        for (level, pop1, pop2) in traverse(self.guide_tree):

            n_pops = len(pop1) + len(pop2)
            haplotypes1 = haplotypes[ np.isin(groups, pop1) ]
            haplotypes2 = haplotypes[ np.isin(groups, pop2) ]

            if len(haplotypes1) < 5:
                cerr('[I - insufficient population size for %s]' % pop1)
            if len(haplotypes2) < 5:
                cerr('[I - insufficient population size for %s]' % pop2)

            # convert haplotypes to allele counts
            ac1 = count_allele(haplotypes1)
            ac2 = count_allele(haplotypes2)

            # calculate highest FST
            FST = []
            num, den = allel.hudson_fst(ac1, ac2)

            # NOTE: the line below might produce warning (invalid value in true_divide)
            # if den == 0, which should be perfectly ok for FST calculation
            fst = num/den

            fst[ np.isnan(fst) ] = 0
            sortidx = np.argsort( fst )

            # get highest FST
            highest_fst_pos = sortidx[-(k+1):-1]
            highest_fst_val = fst[ highest_fst_pos ]
            #cerr('[I - highest FST: %5.4f at %d for pops %s and %s' % (highest_fst_val, highest_fst_pos, pop1, pop2))

            # check suitability of SNPs
            if highest_fst_val.max() < self.min_fst:

                snplist, F = self.select_2(haplotypes1, haplotypes2)
                if snplist:
                    self.log('F: %5.4f SNP: %d for pop %s <> %s' % (F, len(snplist), pop1, pop2))

                    for p in snplist:
                        candidate_L.append( (p, level, n_pops) )
                    continue

                # 2nd approach: find 2 SNPs with highest r^2(st) eg r^2 subpopulation vs r^2 total population
                else:
                    self.log('low FST = %5.4f for %s vs %s' % ( highest_fst_val.max(), pop1, pop2))

            # append to candidate_L
            for p in highest_fst_pos:
                candidate_L.append( (p, level, n_pops) )

        # process candidate_L
        L = np.unique( np.array( sorted( [ x[0] for x in candidate_L] ) ) )

        # return snp position
        return (L, None)


    def select_2(self, haplotypes1, haplotypes2):
        return None, None


class HHFSTDTSelector(HierarchicalFSTSelector):

    code = 'hfst+dt'

    def select_2(self, haplotypes1, haplotypes2):
        """ return (snplist, F):
            snplist - a list of SNP positions after further selection
            F = F score for these particular SNP set """

        X_train =  np.append(haplotypes1, haplotypes2, axis=0)
        y_train =  np.array( [1] * len(haplotypes1) + [2] * len(haplotypes2) )

        best_score = (-1, None, None, None)
        for i in range(3):

            classifier = DecisionTreeClassifier(class_weight='balanced', random_state = self.randomstate)
            classifier = classifier.fit(X_train, y_train)
            features = classifier.tree_.feature

            # remove features with negative position and redundant
            features = np.unique(features[ features >= 0])

            model = FixSNPSelector(features)
            lk_predictions, snplist, _ = fit_and_predict(model, X_train, y_train, X_train, len(features))
            scores = lkprof.calculate_scores(y_train,  lk_predictions, len(features), 'dt', i)

            f_score = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
            if f_score > best_score[0]:
                best_score = (f_score, scores, None, features.tolist())

        return best_score[3], best_score[0]


def prepare_stratified_samples(haplotypes, group_keys, k_fold, haplotype_func=None):
    """ check the suitability of sample sets and modify haplotypes and group_keys properly """

    groups = []
    for group_key, count in zip( *np.unique(group_keys, return_counts=True)):
        # we make sure that every group has at least 2 * k_fold member
        if count < k_fold * 2:
            groups.append( (group_key, math.ceil(k_fold*2 / count)) )

    if len(groups) == 0:
        # nothing to modify
        return (haplotypes, group_keys)

    cerr('[I - prepare_stratified_sample() replicated group: %s]' % ' '.join( x[0] for x in groups ))
    #import IPython; IPython.embed()
    new_haplotypes = [ haplotypes ]
    new_group_keys = [ group_keys ]
    for group_key, m_factor in groups:
        indexes = np.where( group_keys == group_key )
        for i in range(m_factor):
            new_haplotypes.append( haplotypes[indexes] )
            new_group_keys.append( group_keys[indexes] )

    haplotypes = np.concatenate( new_haplotypes, axis=0 )
    group_keys = np.concatenate( new_group_keys, axis=0 )

    return (haplotypes, group_keys)


def fit_and_predict(model, X_train, y_train, X_test, k):
    """ return prediction from likelihood and prediction from original method (if available) as
        (lk_prediction, actual_k, model_prediction)
    """

    # select SNPs based on model
    L, orig_predictions = model.select(X_train, y_train, X_test, k)
    X_train_0 = X_train[:,L]
    X_test_0 = X_test[:,L]

    # create, train and predict using LikelihoodProfile

    lkp = lkprof.LikelihoodProfile()
    lkp.fit(X_train_0, y_train)
    return lkp.predict(X_test_0), L, orig_predictions


def validator_worker( args ):

    """ validator: returns (r, scores, snplist, log)
        where:
            r: repeat identifier
            scores: Panda dataframe containing all scores
            snplist: a dictionary of simid: snplist
            log: list of log message
    """

    model, y, k_list, fold, iteration, simid = args
    pid = os.getpid()

    cerr('[I - pid %d: validator_worker() started]' % pid)

    np.random.seed( simid % pid )
    model.reseed( simid )

    if var_dict['X_shape'] == None:
        X = var_dict['X']
    else:
        cerr('[I - pid %d: validator_worker() is mapping numpy array]' % pid)
        X = np.frombuffer(var_dict['X'], dtype=np.int8).reshape(var_dict['X_shape'])

    # check for sample size suitability for k-folding
    X, y = prepare_stratified_samples( X, y, fold )

    skf = StratifiedKFold(n_splits = fold, shuffle=True, random_state = np.random.randint(1e8))

    results = []
    snps = {}
    k_fold = -1

    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        k_fold += 1

        for k in k_list:

            # best score will be based on highest min F score
            best_score = (-1, None, None, None)
            for i in range(iteration):
                # the iteration here is used for stochastic models where each iteration can yield
                # different result
                lk_predictions, snplist, orig_predictions = fit_and_predict(model, X_train, y_train, X_test, k)
                scores = lkprof.calculate_scores(y_test,  lk_predictions, len(snplist), k, 'lk', simid, k_fold)
                if orig_predictions is not None:
                    orig_scores = lkprof.calculate_scores(y_test, orig_predictions
                            , len(snplist), k, model.code, simid, k_fold)
                else:
                    orig_scores = None

                f_score = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
                if f_score > best_score[0]:
                    best_score = (f_score, scores, orig_scores, snplist.tolist())

            results.append( best_score[1] )
            if best_score[2] is not None:
                results.append( best_score[2] )
            snps['%d/%d/%d/%d' % (simid, k_fold, k, len(best_score[3]))] = best_score[3]

    # reformat model log
    log = [ '[I - {%d} %s]' % (simid, line) for line in model.get_loglines()]

    return (simid, pd.concat( results ), snps, log)


# global variable for multiprocessing
var_dict = {}

def init_worker(X, X_shape):
    var_dict['X'] = X
    var_dict['X_shape'] = X_shape


def validate( model, haplotypes, group_keys, k_list, repeats, fold, iteration
        , outfile, procs=1, outsnp=None, logf=None):

    seed = np.random.randint(1e7)
    arguments = [ (model, group_keys, k_list, fold, iteration, seed+n) for n in range(repeats) ]


    if procs > 1:
        # perform multiprocessing

        # create a shared-memory for numpy array
        cerr('[I - preparing for shared-memory numpy array]')
        X_shape = haplotypes.shape
        X = RawArray('b', X_shape[0]*X_shape[1])
        X_np = np.frombuffer(X, dtype=np.int8).reshape(X_shape)

        # fill share-memory with haplotypes
        np.copyto(X_np, haplotypes)

        with Pool(procs, initializer=init_worker, initargs=(X, X_shape)) as pool:
            c = 0
            for (n, result, snps, log) in pool.imap_unordered(validator_worker, arguments):
                c += 1
                cerr('[I - receiving result from repeat #%d (%d/%d) with %d results]'
                        % (n-seed+1, c, repeats, len(result)))
                # write to temporary files
                if outfile:
                    with open('%s.%d' % (outfile, n), 'wb') as fout:
                        pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)
                if outsnp:
                    with open('%s.%d' % (outsnp, n), 'wb') as fout:
                        pickle.dump(snps, fout, pickle.HIGHEST_PROTOCOL)

                # write to log
                if logf and log:
                    logf.write( '\n'.join( log ) )
                    logf.write( '\n' )

    else:

        init_worker( haplotypes, None )
        c = 0
        for (n, result, snps, log) in map(validator_worker, arguments ):
            c += 1
            cerr('[I - receiving result from repeat #%d (%d/%d) with %d results]'
                    % (n-seed+1, c, repeats, len(result)))

            # write to temporary files
            if outfile:
                with open('%s.%d' % (outfile, n), 'wb') as fout:
                    pickle.dump(result, pickle.HIGHEST_PROTOCOL)
            if outsnp:
                with open('%s.%d' % (outsnp, n), 'wb') as fout:
                    pickle.dump(snps, pickle.HIGHEST_PROTOCOL)

            # write to log
            if logfh and log:
                logfh.write( '\n'.join( log ) )
                logfh.write( '\n' )

    cerr('[I - combining output files]')
    results = []
    snp_tables = {}

    for n in range(repeats):

        filename = '%s.%d' % (outfile, n+seed)
        with open(filename, 'rb') as fin:
            results.append( pickle.load(fin) )
        os.remove(filename)

        filename = '%s.%d' % (outsnp, n+seed)
        with open(filename, 'rb') as fin:
            snp_tables.update( pickle.load(fin) )
        os.remove(filename)


    return pd.concat( results ), snp_tables


def get_model(args):

    if args.method == 'rand':
        return RandomSelector()

    elif args.method == 'dt':
        return DecisionTreeSelector()

    elif args.method == 'hfst':
        return HierarchicalFSTSelector( guide_tree = parse_guide_tree( open(args.guidetree) ) )

    elif args.method == 'hfst+dt':
        return HHFSTDTSelector( guide_tree = parse_guide_tree( open(args.guidetree) ) )

    elif args.method == 'list':
        return FixSNPSelector(snpfile = args.snpfile)

    else:
        cexit('ERR: please provide method')


def simulate(args):

    model = get_model(args)
    cerr('[I - simulation model: %s repeats: %d %d-fold iteration: %d k: %s]'
        % (model.code, args.repeats, args.fold, args.iter, args.k))

    logf = open(args.logfile, 'w') if args.logfile else None

    nalt_parser = naltparser.NAltLineParser( args, datatype='nalt')

    # we need group info
    nalt_parser.parse_grouping()
    group_keys = nalt_parser.group_parser.group_keys

    # remove groups whose samples are less than 3
    suitable_groups = set( [ g for g in nalt_parser.group_parser.groups
                                if len(nalt_parser.group_parser.groups[g]) > 2])

    cerr('[I - reading %s samples]' % len(nalt_parser.samples))
    mask = [False] * len(nalt_parser.samples)
    for i in range(len(nalt_parser.samples)):
        if group_keys[i] in suitable_groups:
            mask[i] = True
    samples = list(itertools.compress(nalt_parser.samples, mask))
    group_keys = np.array(list(itertools.compress(group_keys, mask)))

    cerr('[I - masking to %d samples]' % len(samples))
    region = nalt_parser.parse_whole(mask=mask)

    if args.mac > 0:
        cerr('[I - filtering for MAC = %d]' % args.mac)
        region.filter_mac(args.mac, inplace=True)

    cerr('[I - preparing haplotypes]')
    haplotypes = region.haplotypes()

    k_list = [ int(x) for x in args.k.split(',') ]
    results, snp_table = validate( model, haplotypes, group_keys, k_list = k_list
        , repeats = args.repeats, fold = args.fold, iteration = args.iter, outfile=args.outfile
        , procs = args.j, outsnp = args.outsnp, logf = logf )

    results.to_csv(args.outfile, sep='\t', index=False)
    cerr('[I - writing scores to %s]' % args.outfile)

    if args.outsnp:
        import yaml
        yaml.dump(snp_table, open(args.outsnp, 'w'))
        cerr('[I - writing SNP table to %s]' % args.outsnp )

    if logf:
        cerr('[I - writing log to %s]' % args.logfile)


def train(args):

    nalt_parser = naltparser.NAltLineParser( args, dataype='nalt')

    # we need group info
    nalt_parser.parse_grouping()
    group_keys = nalt_parser.group_parser.group_keys


    cerr('[I - masking to %d samples]' % len(samples))
    region = nalt_parser.parse_whole(mask=mask)
    haplotypes = region.haplotypes()

    matrix_profiles = cross_train( haplotypes, group_keys, repeats, fold, )

def predict(args):

    # read input profiles as Python pickle

    profiles = None



mode = {    'simulate': simulate,
            'train': train,
            'predict': predict }

def lkest(args, config=None):

    # check for config as dictionary
    if config:

        # general
        args.mode = config.get('mode', None) or args.mode


    try:
        func = mode[args.mode]
    except KeyError:
        cexit('[E - mode %s does not exist!]' % args.mode)

    # run mode
    start_time = time.monotonic()
    func(args)
    cerr('[I - finished in %6.2f minute(s) at %s]'
            % (time.monotonic() - start_time, datetime.datetime.now()))



def main(args):

    # check if using YAML config file
    config = yaml.load(open(args.config)) if args.config else None

    lkest(args, config)
