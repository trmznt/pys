#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an FST file with the following format:
# AVG_FST MAX_FST POP1 POP2 ...


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import itertools
import allel

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('--grp1', required=True)
    p.add_argument('--grp2', required=True)

    return p


def main( args ):

    geno2pairfst( args )


def geno2pairfst( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.diploid_translator)

    cout('Grouping:')
    groups = lineparser.parse_grouping()
    for k in groups:
        cout(' %12s %3d' % (k, len(groups[k])))

    FST = [] # FST indexed by group_keys
    group_keys = sorted(groups.keys())
    cout(group_keys)

    # gathering groups
    grp1 = list(itertools.chain.from_iterable(
        groups[k] for k in args.grp1.split(',')
    ))
    grp2 = list(itertools.chain.from_iterable(
        groups[k] for k in args.grp2.split(',')
    ))
 
    # output to file
    FST = []


    idx = 0
    for (posinfo, genolist) in lineparser.parse():

        idx += 1
        genoarray = allel.GenotypeArray( [genolist]  )


        # calculate FST per group against other samples

        ac_g1 = genoarray.count_alleles(subpop = grp1)
        ac_g2 = genoarray.count_alleles(subpop = grp2)
        num, den = allel.stats.hudson_fst(ac_g1, ac_g2)
        fst = num[0]/den[0]
        if not (0.0 <= fst <= 1.0):
                fst = 0
        FST.append( (fst, posinfo) )

    FST.sort(reverse=True)
    for fst, posinfo in FST[:10]:
        cout('%s\t%s\t%s\t%5.4f' % (posinfo[0], posinfo[1], posinfo[4], fst))