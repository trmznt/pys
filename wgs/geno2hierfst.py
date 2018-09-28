#!/usr/bin/env spcli

# this command runs hierarchical FST comparison

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import itertools
import allel

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('--hierfile', required=True, help="file describing hierarchical groups")

    return p


def main( args ):

    geno2hierfst( args )


def geno2hierfst( args ):

    genoparser = tabparser.GenotypeLineParser( args )
    genoparser.set_translator(genoparser.diploid_translator)

    cerr('Grouping:')
    groups = genoparser.parse_grouping()
    for k in groups:
        cout(' %12s %3d' % (k, len(groups[k])))

    hierarchy = []
    with open(args.hierfile) as hierfile:
        for line in hierfile:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            partitions = line.split('\t')
            print(partitions)
            par1 = list(itertools.chain.from_iterable(
                    groups[k] for k in partitions[0].split(',')))
            par2 = list(itertools.chain.from_iterable(
                    groups[k] for k in partitions[1].split(',')))
            hierarchy.append( (par1, par2) )
    cerr('[I: preparing %d hierarchy]' % len(hierarchy))

    cerr('[I: reading genotype file...]')
    genotypes = genoparser.parse_all()
    genoarray = allel.GenotypeArray( genotypes )
    #import IPython; IPython.embed()
    del genotypes

    selected_positions = []
    c = 1
    for (grp1, grp2) in hierarchy:
        cerr('[I: processing hierarchy #%d]' % c)
        FST = []
        ac_g1 = genoarray.count_alleles(subpop = grp1)
        ac_g2 = genoarray.count_alleles(subpop = grp2)
        #import IPython; IPython.embed()
        num, den = allel.stats.hudson_fst(ac_g1, ac_g2)
        fst = num/den

        for p,v in zip(genoparser.position, fst):
            if not (0.0 <= v <= 1.0):
                v = 0
            FST.append( (v,p) )
        FST.sort( reverse=True )
        cumulative_fst = 0.0
        for (v, p) in FST:
            if cumulative_fst > 2.0:
                break
            selected_positions.append( (p, v) )
            cumulative_fst += v
        c += 1

    for (p, v) in selected_positions:
        cout('%s\t%s\t%s\t%5.4f' % (p[0], p[1], p[4], v))



