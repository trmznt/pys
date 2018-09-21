#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an FST file with the following format:
# CHROM POS GENE SNPS HAPLOTYPES FST AVG_FST MAX_FST MIN_FST POP1 POP2 ...


from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np
import itertools

from seqpy.core.funcs import pruner

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-s', '--scheme', type=int, default=2)
    p.add_argument('--threshold', type=float, default=0.25)
    p.add_argument('--region', default='chrom')
    p.add_argument('-o', '--outfile', default='outfile.pruned.txt')

    return p


def main( args ):

    geno2prune( args )


def geno2prune( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)

    cerr('I: start parsing genotype file...')

    outfile = open(args.outfile, 'w')
    outfile.write(lineparser.get_position_header())
    outfile.write('\n')

    iter_func = {
        'chrom': lineparser.parse_chromosomes,
        'whole': lineparser.parse_whole,
        'gene': lineparser.parse_genes
    }[args.region]

    for region in iter_func():
        cerr('I: generating genotype array for %s' % region.name)
        genoarray = allel.GenotypeArray( region.genotypes() )

        cerr('I: pruning for %s' % region.name)
        if args.scheme == 1:
            index = pruner.prune_1(genoarray, args.threshold)
        elif args.scheme == 2:
            index = pruner.prune_2(genoarray, region.positions(), args.threshold)
        elif args.scheme == 3:
            index = pruner.prune_3(genoarray, region.positions(), args.threshold)
        else:
            cexit('E: scheme type undefined!')

        new_positions = itertools.compress(region.positions(), index)

        for pos in new_positions:
            outfile.write('%s\n' % '\t'.join( pos ))
