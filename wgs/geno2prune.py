#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an FST file with the following format:
# CHROM POS GENE SNPS HAPLOTYPES FST AVG_FST MAX_FST MIN_FST POP1 POP2 ...


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np

from seqpy.core.funcs import pruner

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-s', '--scheme', type=int, default=2)
    p.add_argument('-o', '--outfile', default='outfile.pruned.txt')

    return p


def main( args ):

    geno2prune( args )


def geno2prune( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.haploid_translator)

    cerr('I: start parsing genotype file...')

    for region in lineparser.parse_chromosomes():
        cerr('I: pruning %s' % region.name)
        genoarray = allel.GenotypeArray( region.genotypes() )
        index = pruner.prune_2(genoarray, region.positions())