# ralt2qc.py
#
# report sample and snp QC



from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np

from seqpy.core.cfuncs import genoutils
from seqpy.core.bioio import naltparser

import random

def init_argparser():
    p = arg_parser("Assess SNP and sample quality from ralt file")


    p = naltparser.init_argparser(p)

    return p

def main(args):
    nalt2qc( args )


def nalt2qc( args ):
    """ write to out.imiss & out.lmiss
        for each sample and SNPS, evaluate:
            N_MISS
            F_MISS
            N_HETS
            F_HETS

        out.imiss:
        SAMPLE N_SNP N_MISS F_MISS N_HETS F_HETS

        out.lmiss:
        CHR POS N_SNP N_MISS F_MISS N_HETS F_HETS


     """

    nalt_parser = naltparser.NAltLineParser( args, datatype='nalt')

    samples = nalt_parser.parse_samples()
    whole = nalt_parser.parse_whole()

    # create an array for N samples with column:
    # N_SNP N_MISS N_HETS

    # for SNPs, appending this list
    snps = []

    for (pos, n_alt) in whole.parse_positions():
        # do the counting here

    # do the stats here

    # write the outputs






