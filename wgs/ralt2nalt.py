# ralt2nalt.py
#
# convert ralt file to nalt file


from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

import numpy as np

from seqpy.core.cfuncs import genoutils
from seqpy.core.bioio import naltparser

def init_argparser():
    p = arg_parser("Convert r_alt datafile to n_alt datafile")
    p.add_argument('-o', '--outfile', default='outdata')
    p.add_argument('-r', '--hetratio', type=float, default=0.25,
                    help = "ratio for calling heterozygosity")
    p.add_argument('-m', '--major', action = "store_true",
                    help = "set all caling to major allele")
    p.add_argument('--missingtohet', action = 'store_true',
                    help = "set missing as heterozygosity")
    p.add_argument('--outfmt', default='npy', choices=['npy', 'pickle', 'tab'])
    p.add_argument('--autofilename', default=True, action='store_true')
    # option to treat all missing as hets

    p = naltparser.init_argparser(p, with_group=False, with_position=False)

    return p

def main(args):
    ralt2nalt( args )


def ralt2nalt( args ):

    ralt_parser = naltparser.NAltLineParser( args, datatype='ralt'
        , with_group=False, with_position=False)

    region = ralt_parser.parse_whole()

    # convert to n_alt
    cerr('[I - converting to nalt format]')
    cerr( '[ M dtype: {}]'.format(region.M.dtype) )
    region.ralt_to_nalt(hetratio = args.hetratio if not args.major else -1)
    cerr('[ M dtype: {}]'.format(region.M.dtype) )

    region.save(args.outfmt, prefixname=args.outfile, autofilename=args.autofilename
            , with_position=False)
    return

    # write to outfile
    with open(args.outfile, 'w') as outfile:
        # write header
        outfile.write(ralt_parser.get_sample_header())
        outfile.write('\n')
        np.savetxt(outfile, region.M, fmt ='%d', delimiter = '\t')

    cerr('[I: finish writing to %s' % args.outfile)

