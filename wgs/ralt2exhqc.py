#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser

import numpy as np
import pandas as pd

import io, time


def init_argparser(p=None):
    if not p:
        p = arg_parser("ralt2allqc - checking exhaustively")

    p = naltparser.init_argparser(p, with_group=False, with_position=False)

    p.add_argument('--mac', type=int, default=1, help="Minor Allele Count")
    p.add_argument('--miss', type=float, default=1.0, help='Missingness fraction threshold')
    p.add_argument('--mindepth', type=int, default=0, help='minimum depth, need to use depth file')
    p.add_argument('-k', type=int, default=-1, help='start from individual/samples no K')
    p.add_argument('-o', '--outfile', default='out.exhqc.txt')
    p.add_argument('-s', type=int, default=-1, help='no of sample to be exported')

    return p


def main( args ):

    ralt2exhqc( args )


def ralt2exhqc( args ):

    cerr('[I - reading input files]')

    start_time = time.monotonic()
    M, sample_idx, site_idx = read_data( args )
    cerr('[I - reading %d sites for %d samples in %d secs]'
        % (len(site_idx), len(sample_idx), time.monotonic() - start_time))

    # create sample missingness
    indv_missing = np.count_nonzero(M < args.mindepth, axis=0) / len(sample_idx)
    complete_indv = np.count_nonzero(indv_missing == 0.0)

    # sorted by args
    indv_idx = np.argsort( indv_missing )

    if args.s > 0:
        do_export_ralt(M, sample_idx, site_idx, indv_idx, args)
    else:
        do_exhqc_ralt(M, sample_idx, site_idx, indv_idx, complete_indv, args)


def do_export_ralt(M, sample_idx, site_idx, indv_idx, args):

    cerr('[I - exporting sample and position indexes for %d samples' % args.s)
    collected_samples = indv_idx[:args.s]
    filt_site_idx, inf_site_idx = filter_site_idx(M, collected_samples, site_idx, mac=args.mac)
    np.savetxt('exhqc.indv.txt', collected_samples, fmt='%d')
    np.savetxt('exhqc.pos.txt', inf_site_idx, fmt='%d')


def do_exhqc_ralt(M, sample_idx, site_idx, indv_idx, complete_indv, args):

    n_samples = []
    n_snps = []
    n_isnps = []

    start_indv = len(indv_idx) if args.k < 0 else args.k
    cerr('[I - iteration start from %d to %d samples' % (start_indv, complete_indv-1))
    for k in range(start_indv, complete_indv-1, -1):
        sample_indexes = indv_idx[:k]

        cerr('[I - filtering using %d samples]' % len(sample_indexes))
        filt_site_idx, inf_site_idx = filter_site_idx(M, sample_indexes, site_idx, args.miss, args.mac)
        n_samples.append( len(sample_indexes) )
        n_snps.append( len(filt_site_idx) )
        n_isnps.append( len(inf_site_idx))

    df = pd.DataFrame({ 'N_SAMPLES': n_samples, 'N_SNPS': n_snps, 'N_ISNPS': n_isnps})
    df.to_csv(args.outfile, sep='\t', index=False)


def filter_site_idx(M, sample_idx, site_idx, lmiss=1.0, mac=1):
    """  return new site_idx where sites are complete with proper MAC """

    M = M[:, sample_idx]

    M2, site_idx1, _ = filter_lmiss(M, site_idx, sample_idx, lmiss)
    if len(site_idx1) == 0:
        return site_idx1, site_idx1

    M2, site_idx2, _ = filter_mac(M2, site_idx1, sample_idx, mac)

    return site_idx1, site_idx2


def check_sanity(M, site_idx, sample_idx):
    shape = M.shape
    # sanity checking
    if len(sample_idx) != shape[1]:
        print( len(sample_idx), shape[1] )
        cexit('[E - inconsistent M shape and no of samples!]')
    if len(site_idx) != shape[0]:
        print( len(site_idx), shape[0] )
        cexit('[E - inconsistent M shape and no of sites!]')


def filter_lmiss(M, site_idx, sample_idx, lmiss):

    cerr('[I - filtering for SNP missingness < %4.3f]' % lmiss)
    check_sanity(M, site_idx, sample_idx)
    site_missingness = np.count_nonzero(M < 0, axis=1) / len(sample_idx)
    indexes = np.where( site_missingness <= (1.0 - lmiss) )

    #if len(indexes[0]) > 0:
    #    import IPython; IPython.embed()

    M2 = M[ indexes[0], : ]
    site_idx2 = site_idx[ indexes[0] ]
    cerr('[I - missingness filter kept  %d from %d sites]' % (len(site_idx2), len(site_idx)))
    #import IPython; IPython.embed()

    return M2, site_idx2, sample_idx


def filter_mac(M, site_idx, sample_idx, mac):

    cerr('[I - filtering for MAC >= %d]' % mac)
    check_sanity(M, site_idx, sample_idx)
    allele_0 = np.count_nonzero(M < 0.5, axis=1)
    allele_1 = len(sample_idx) - allele_0
    allele_mac = np.minimum(allele_0, allele_1)
    indexes = np.where( allele_mac >= mac )

    M2 = M[indexes[0], :]
    site_idx2 = site_idx[ indexes[0] ]
    cerr('[I - MAC filter kept %d from %d sites]' % (len(site_idx2), len(site_idx)))
    #import IPython; IPython.embed()

    return M2, site_idx2, sample_idx


def read_data_xxx( args ):
    """ return M, sample_idx, site_idx """
    df = pd.read_csv(args.infile, sep='\t', dtype=float,
            nrows=args.n if args.n > 0 else None)
    samples = df.columns
    sample_idx = np.arange(len(samples))
    M = np.rint(df.values).astype(np.short)
    site_idx = np.arange(len(M))

    return M, sample_idx, site_idx


def read_data( args ):

    nalt_parser = naltparser.NAltLineParser( args, with_group=False, with_position=False)
    region = nalt_parser.parse_whole()
    return region.M, np.arange(len(region.M[0])), np.arange(len(region.M))

