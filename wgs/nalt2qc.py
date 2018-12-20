""" ralt2qc.py
    report sample and SNP QC
"""
import random, time

import numpy as np

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import genoutils


def init_argparser():
    parser = arg_parser('Assess SNP and sample quality from nalt file')
    parser = naltparser.init_argparser(parser)

    return parser


def nalt2qc(args):
    """ write to out.imiss & out.lmiss
        for each sample and SNPS, evaluate:
            N_MISS
            F_MISS
            N_HETS
            F_HETS

        out.imiss:
        SAMPLE N_SNP N_MISS F_MISS N_HETS F_HETS

        out.lmiss:
        CHR POS N_SAMPLE N_MISS F_MISS N_HETS F_HETS
    """

    start_time = time.monotonic()
    nalt_parser = naltparser.NAltLineParser(args, datatype='nalt')
    samples = nalt_parser.parse_samples()
    whole = nalt_parser.parse_whole()
    cerr('[I - reading input file in %s secs]' % ( time.monotonic() - start_time))

    # create an array for N samples with column:
    # N_SNP N_MISS N_HETS
    asamples = np.zeros(shape=(len(samples), 3))

    # container for lmiss output
    chr_pos = []
    snps = []

    # -1 indicate missing data, 1 indicate heterozygous SNP
    for pos, n_alt in whole.parse_positions():
        # gather imiss data
        asamples[range(n_alt.size), 0] += 1
        asamples[n_alt == -1, 1] += 1
        asamples[n_alt == 1, 2] += 1

        # gather lmiss data
        chromosome = pos[0]
        position = pos[1]
        n_sample = len(n_alt)
        n_miss = np.where(n_alt == -1)[0].size
        n_het = np.where(n_alt == 1)[0].size
        chr_pos.append((chromosome, position))
        snps.append((n_sample, n_miss, n_het))

    # create imiss stats
    imiss = np.zeros(shape=(len(samples), 5))
    n_snps = asamples[:, 0].sum()
    imiss[:, 0] = asamples[:, 0]
    imiss[:, 1] = asamples[:, 1]
    imiss[:, 2] = asamples[:, 1] / n_snps
    imiss[:, 3] = asamples[:, 2]
    imiss[:, 4] = asamples[:, 2] / n_snps

    # create lmiss stats
    snps = np.array(snps)
    n_samples = snps[:, 0].sum()
    lmiss = np.zeros(shape=(len(snps), 5))
    lmiss[:, 0] = snps[:, 0]
    lmiss[:, 1] = snps[:, 1]
    lmiss[:, 2] = snps[:, 1] / n_samples
    lmiss[:, 3] = snps[:, 2]
    lmiss[:, 4] = snps[:, 2] / n_samples

    # output imiss
    with open('out.imiss', 'w') as iout:
        iout.write('INDV\tN_SNP\tN_MISS\tF_MISS\tN_HETS\tF_HETS\n')
        for sample, i in zip(samples, imiss):
            iout.write('{}\t'.format(sample))
            iout.write('{}\n'.format('\t'.join(map(str, i))))

    with open('out.lmiss', 'w') as lout:
        lout.write('CHROM\tPOS\tN_SAMPLE\tN_MISS\tF_MISS\tN_HETS\tF_HETS\n')
        for (chrom, pos), l in zip(chr_pos, lmiss):
            lout.write('{}\t{}\t'.format(chrom, pos))
            lout.write('{}\n'.format('\t'.join(map(str, l))))


# the arguments are parsed by spcli if there is a init_argparser
def main(args):
    nalt2qc(args)
