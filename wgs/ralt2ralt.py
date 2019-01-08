#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser

import numpy as np

import io, time


def init_argparser(p=None):
    if not p:
        p = arg_parser("ralt2allqc - checking exhaustively")

    p = naltparser.init_argparser(p)

    p.add_argument('--indvindex', default=None)
    p.add_argument('--posindex', default=None)
    p.add_argument('--type', default='ralt')
    p.add_argument('-o', '--outfile', default='outfile')

    return p


def main( args ):

    ralt2ralt( args )


def read_data( args ):
    """ return M, sample_idx, site_idx """
    df = pd.read_csv(args.infile, sep='\t',
            dtype=float if args.type=='ralt' else int)
    samples = df.columns
    sample_idx = np.arange(len(samples))
    M = np.rint(df.values).astype(np.short)
    site_idx = np.arange(len(M))

    return M, sample_idx, site_idx



def ralt2ralt( args ):

    cerr('[I - reading input files]')

    start_time = time.monotonic()

    alt_parser = naltparser.NAltLineParser(args, datatype=args.type)

    whole_region = alt_parser.parse_whole()
    samples = alt_parser.parse_samples()

    if args.posindex:
        pos_indexes = np.loadtxt(args.posindex, dtype=int)
        whole_region.filter_positions(pos_indexes)

    if args.indvindex:
        indv_indexes = np.loadtxt(args.indvindex, dtype=int)
        whole_region.filter_samples(indv_indexes)
        samples = samples[indv_indexes]

    # save to outfile

    outmatrix = args.outfile + ('.ralt.txt' if args.type == 'ralt' else '.nalt.txt')
    outpos = args.outfile + '.pos.txt'

    whole_region.df_M.to_csv(outmatrix, sep='\t', index=False)
    whole_region.df_P.to_csv(outpos, sep='\t', index=False)
    cerr('[I - writing to file: %s and %s' % (outmatrix, outpos))

    return

    with open(outmatrix, 'wt') as f_matrix, open(outpos, 'wt') as f_pos:
        f_matrix.write('\t'.join(samples))
        f_matrix.write('\n')
        np.savetxt(f_matrix, whole_region.M, delimiter='\t',
                fmt='%4.3f' if args.type == 'ralt' else '%d')

        f_pos.write('\t'.join(alt_parser.position_parser.header))
        np.savetxt(f_pos, whole_region.P, delimiter='\t')

    cerr('[I - writing to file: %s and %s' % (outmatrix, outpos))
