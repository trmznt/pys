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

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.fstgene.txt')

    return p


def main( args ):

    geno2fstgene( args )


def geno2fstgene( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.diploid_translator)

