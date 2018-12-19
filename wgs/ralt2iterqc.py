#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import io, time


def init_argparser(p=None):
    if not p:
        p = arg_parser("ralt2iterqc - iter")
    p.add_argument('--filter', default='')
    p.add_argument('--lmiss', type=float, default=-1)
    p.add_argument('--imiss', type=float, default=-1)
    p.add_argument('--mac', type=int, default=-1)
    p.add_argument('--iter', type=int, default=1)
    p.add_argument('-n', type=int, default=-1)
    p.add_argument('infile')

    return p


def main( args ):

    ralt2iterqc( args )


def ralt2iterqc( args ):

    cerr('[I - reading input files]')

    start_time = time.monotonic()
    df = pd.read_csv(args.infile, sep='\t', dtype=float,
            nrows=args.n if args.n > 0 else None)
    samples = df.columns
    sample_idx = np.arange(len(samples))
    M = df.values
    site_idx = np.arange(len(M))

    cerr('[I - reading %d sites for %d samples in %d secs]'
        % (len(site_idx), len(sample_idx), time.monotonic() - start_time))

    for i in range(args.iter):
        cerr('[I - ITER -> %d]' % (i+1))
        if args.lmiss > 0:
            M, site_idx, sample_idx = filter_lmiss(M, site_idx, sample_idx, args.lmiss)
        if args.imiss > 0:
            M, site_idx, sample_idx = filter_imiss(M, site_idx, sample_idx, args.imiss)
        if args.mac > 0:
            M, site_idx, sample_idx = filter_mac(M, site_idx, sample_idx, args.mac)


    # filtering order: lmiss > imiss > mac


def check_sanity(M, site_idx, sample_idx):
    shape = M.shape
    # sanity checking
    if len(sample_idx) != shape[1]:
        print( len(sample_idx), shape[1] )
        cexit('[E - inconsistent M shape and no of samples!]')
    if len(site_idx) != shape[0]:
        cexit('[E - inconsistent M shape and no of sites!]')


def filter_lmiss(M, site_idx, sample_idx, lmiss):

    cerr('[I - filtering for SNP missingness < %4.3f]' % lmiss)
    check_sanity(M, site_idx, sample_idx)
    site_missingness = np.count_nonzero(M < 0, axis=1) / len(sample_idx)
    indexes = np.where( site_missingness < (1.0 - lmiss) )

    M2 = M[ indexes ]
    site_idx2 = site_idx[ indexes ]
    cerr('[I - keeping %d from %d sites]' % (len(site_idx2), len(site_idx)))
    #import IPython; IPython.embed()

    return M2, site_idx2, sample_idx


def filter_imiss(M, site_idx, sample_idx, imiss):

    cerr('[I - filtering for sample missingness < %4.3f]' % imiss)
    check_sanity(M, site_idx, sample_idx)
    indv_missingness = np.count_nonzero(M < 0, axis=0) / len(site_idx)
    indexes = wp.where( indv_missingness < (1.0 - imiss) )

    M2 = M[:, indexes]
    sample_idx2 = sample_idx[ indexes ]
    cerr('[I - keeping %d from %d samples]' % (len(sample_idx2), len(sample_idx)))

    return M, site_idx, sample_idx


def filter_mac(M, site_idx, sample_idx, mac):

    cerr('[I - filtering for MAC > %d]' % mac)
    check_sanity(M, site_idx, sample_idx)
    allele_0 = np.count_nonzero(M < 0.5, axis=1)

    return M, site_idx, sample_idx
