#
# given sample pw dist & quality list, provide a list of samples to be removed

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser

import numpy as np, pandas as pd

from types import SimpleNamespace

def init_argparser(p=None):

    if not p:
        p = arg_parser('dist2clonalqc.py')
    p.add_argument('--qualfile', default=None)
    p.add_argument('--picklefile', default=None)
    p.add_argument('--threshold', type=float, default=0.01/100)
    p.add_argument('-o', '--outfile', default=None)
    p.add_argument('infile')
    return p


def clonal_index(dist, qual, samples, threshold=0.01/100):

    clonals = []
    for i,j in zip( *np.where( dist <= threshold ) ):
        if i == j: continue
        if i in clonals or j in clonals:
            continue
        cerr('[I - clonal %s <> %s : %f]' % (samples[i], samples[j], dist[i,j]))
        if qual[i] < qual[j]:
            clonals.append(i)
        else:
            clonals.append(j)

    return clonals

def dist2clonalqc( args ):

    # read distance matrix
    df = pd.read_csv(args.infile, sep='\t')
    samples = df.columns
    D = df.values

    # read quality file or pickled ralt/nalt file
    if args.picklefile:

        nalt_args = SimpleNamespace( infile = args.picklefile, fmt = 'pickle', n = -1)
        nalt_parser = naltparser.NAltLineParser(nalt_args
                , with_group=False, with_position=False)
        region = nalt_parser.parse_whole()
        qual = np.count_nonzero( region.M == -1, axis = 0 )

    else:
        cexit('ERR: other input file has not been defined')

    clonal_samples = clonal_index(D, qual.max() - qual, samples, args.threshold)
    cerr('[I - removing %d clonal samples]' % len(clonal_samples))
    if args.outfile:
        np.savetxt(args.outfile, samples[clonal_samples], fmt='%s')


def main( args ):

    dist2clonalqc( args )