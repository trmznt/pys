# ralt2nalt.py
#
# convert ralt file to nalt file


from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np

from seqpy.core.cfuncs import genoutils
from seqpy.core.bioio import naltparser

import random

def init_argparser():
    p = arg_parser("Evaluate a set of SNPs for its likelihood-based assignment")
    p.add_argument('-k', default='1000,500,100', help='k list')

    p = naltparser.init_argparser(p)

    return p

def main(args):
    nalt2lktest( args )


def nalt2lktest( args ):

    nalt_parser = naltparser.NAltLineParser( args, datatype='nalt')

    # we need group info
    nalt_parser.parse_grouping()

    region = nalt_parser.parse_whole()

    overall_accuracies = []
    for k in [ int(x) for x in args.k.split(',')]:

        cerr('[I: random sites for k=%d]' % k)
        accuracies = []
        for i in range(100):
            rand_reg = random_sites(region, k)
            tester = LikelihoodTester(rand_reg, nalt_parser.samples, nalt_parser.group_parser)
            tester.prepare_matrix()
            accuracies.append(tester.evaluate())
            cerr('[I: iteration of %d]' % i)

        overall_accuracies.append( (k, accuracies) )

    import IPython; IPython.embed()


class LikelihoodTester(object):

    def __init__(self, region, samples, groups):
        self.region = region
        self.groups = groups
        self.samples = samples
        self.matrices = None


    def estimate(self, n_alt_data):

        lks = self.calculate_likelihoods(n_alt_data)
        lks.sort()
        return lks[-1][1]

    def evaluate(self):
        """ return the accuracy of the region """

        l = len(self.region.M)
        match = mismatch = 0
        for i in range(len(self.samples)):
            # create sample array
            n_alt_data = np.zeros(l)
            for x in range(l):
                n_alt_data[x] = self.region.M[x][i]

            est_group = self.estimate(n_alt_data)

            if est_group == self.groups.group_info[self.samples[i]]:
                match += 1
            else:
                mismatch += 1

        return (match / (match+mismatch))


    def prepare_matrix(self):
        """ create a set of matrices for each group """

        # create empty matrices
        M = {}
        region_M = np.array(self.region.M)
        l = len(region_M)
        for group in self.groups.groups:
            M[group] = np.full(shape=(l, 2), fill_value=0.01)

        # filling-up matrices
        # we are iterating column wise
        for i, sample in enumerate(self.samples):
            group = self.groups.group_info[sample]
            n_alt = region_M[:, i]
            M_group = M[group]

            M_group[n_alt == 0, 0] += 2
            M_group[n_alt == 1] += 1
            M_group[n_alt == 2, 1] += 2

        # convert M into fractions
        for group in M:
            M_group = M[group]
            total = np.sum(M_group, axis=1)
            M[group] = np.log(M_group / total[:, None])

        self.matrices = M
        return M


    def calculate_likelihoods(self, n_alt_data):
        """ calculate log likelihood of each group """

        R = []
        for g in self.matrices:
            R.append(
                ( self.calculate_likelihood(self.matrices[g], n_alt_data), g)
            )

        return R


    def calculate_likelihood(self, M, n_alt_data):
        n_alt_0 = 2 * M[np.where(n_alt_data == 0), 0]
        n_alt_1 = M[np.where(n_alt_data == 1), 0] + M[np.where(n_alt_data == 1), 1]
        n_alt_2 = 2 * M[np.where(n_alt_data == 2), 1]
        log_lk = np.sum(n_alt_0) + np.sum(n_alt_1) + np.sum(n_alt_2)

        return log_lk


def random_sites(region, L):
    """ creates a new region with randomly selected positions """

    new_region = naltparser.Region('random')
    for i in range(L):
        new_region.M.append( random.choice(region.M))
    return new_region
