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
    p.add_argument('--includepos', default=None)
    p.add_argument('--excludesample', default=None)
    p.add_argument('--includesample', default=None)
    p.add_argument('--mac', default=0, type=int)
    p.add_argument('--type', default='ralt')
    p.add_argument('--autofilename', default=False, action='store_true')
    p.add_argument('--outfmt', default='text', choices=['text', 'pickle', 'npy'])
    p.add_argument('-o', '--outfile', default='outfile')

    return p


def main( args ):

    ralt2ralt( args )


def read_data( args ):
    """ Note: this function is no longer used;
        return M, sample_idx, site_idx
    """
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
        cerr('[I - filtering for %d SNP position]' % len(pos_indexes))
        whole_region.filter_positions(pos_indexes)
    elif args.includepos:
        with open(args.includepos) as f_posline:
            poslines = [ x.split() for x in f_posline ]
            if poslines[0][0] == 'CHROM' and poslines[0][1] == 'POS':
                del poslines[0]
        whole_region.filter_poslines(poslines, inplace=True)

    if args.indvindex:
        indv_indexes = np.loadtxt(args.indvindex, dtype=int)
        whole_region.filter_samples(indv_indexes)
        samples = samples[indv_indexes]

    if args.excludesample:
        excluded_samples = np.loadtxt(args.excludesample, dtype=str)
        excluded_indexes = np.where(np.array(samples) == excluded_samples[:,None])[1]
        indv_indexes = np.array( list(set( range(len(samples))) - set(excluded_indexes)) )
        cerr('[I - excluding %d samples]' % len(excluded_indexes))
        whole_region.filter_samples(indv_indexes)
        samples = samples[indv_indexes]

    if args.includesample:
        included_samples = np.loadtxt(args.includesample, dtype=str)
        included_indexes = np.where(np.array(samples) == included_samples[:,None])[1]
        cerr('[I - including {} | {} out of {} samples]'.format(
            len(included_samples), len(included_indexes), len(samples)))
        whole_region.filter_samples(included_indexes)
        samples = samples[included_indexes]

    if args.mac > 0:
        whole_region.filter_mac(args.mac)

    # save to outfile
    whole_region.save(args.outfmt, prefixname=args.outfile, autofilename=args.autofilename
            , with_position=True)
    return


    if args.autofilename:
        args.outfile = '%s-%d-%d' % (
                'r' if args.type == 'ralt' else 'n',
                len(samples), len(whole_region.M)
        )

    cerr('[I - writing to outfiles]')
    outmatrix = args.outfile + ('.ralt' if args.type == 'ralt' else '.nalt')
    if args.outfmt == 'pickle':
        outmatrix = outmatrix + '.pickle.gz'
        whole_region.df_M.to_pickle(outmatrix)
    else:
        outmatrix = outmatrix + '.txt.gz'
        whole_region.df_M.to_csv(outmatrix, sep='\t', index=False)

    outpos = args.outfile + '.pos.txt.gz'
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
