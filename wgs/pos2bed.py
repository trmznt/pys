#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

# the BED numbering scheme is based from this URL:
# https://www.biostars.org/p/84686/
#
# diagram:
# chr1  T A C G T
#      |||||||||||
# 1-   |1|2|3|4|5|
# 0-   0 1 2 3 4 5
#
# 1-based: VCF, SAM, GFF, posfile
# 0-based: BED

def init_argparser(p=None):

    if not p:
        p = arg_parser("pos2bed - convert SNP position file to BED format")

    p.add_argument('-o', '--outfile', default='out.pos2bed.txt')
    p.add_argument('infile')

    return p

def main( args ):

    pos2bed( args )


def pos2bed( args ):


    # read posfile as simple tab-delimited file
    positions = []
    with open(args.infile) as fin:
        next(fin)
        for line in fin:
            t = line.split('\t')
            positions.append( ( t[0], int(t[1]), t[2]) )

    # convert to BED

    with open(args.outfile, 'w') as fout:
        for seq, pos, ref in positions:
            fout.write( '%s\t%d\t%d\t%s_%d_%s\n' % (seq, pos-1, pos, seq, pos, ref))

    cerr('[I - writing output to %s]' % args.outfile )