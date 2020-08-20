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

    p.add_argument('--namecol', type=int, default=-1)
    p.add_argument('--refcol', type=int, default=-1)
    p.add_argument('--microhaps', action='store_true', default=False)
    p.add_argument('-o', '--outfile', default='out.pos2bed.txt')
    p.add_argument('infile')

    return p

def main( args ):

    pos2bed( args )


def pos2bed( args ):


    # read posfile as simple tab-delimited file
    positions = []
    with open(args.infile) as fin:
        header = next(fin)
        for line in fin:
            positions.append( line.strip().split('\t') )

    # check options
    if args.namecol < 0 and args.refcol < 0:
        # use standar format CHROM POS REF
        args.refcol = 2

    if args.microhaps:
        return pos2bed_microhaps(args, positions)

    # convert to BED
    with open(args.outfile, 'w') as fout:
        last_name = ''
        last_idx = 0
        for entry in positions:
            seq = entry[0]
            pos = int(entry[1])
            ref = entry[args.refcol] if args.refcol > 0 else 'N'
            name = entry[args.namecol] if args.namecol > 0 else ('%s_%d_%s' % (seq, pos, ref))
            if name == last_name:
                last_idx += 1
                name = '%s_%d' % (name, last_idx)
            else:
                last_name = name
                last_idx = 0
            fout.write( '%s\t%d\t%d\t%s\n' % (seq, pos-1, pos, name))

    cerr('[I - writing SNP-based BED to %s]' % args.outfile )

def pos2bed_microhaps(args, positions):

    if args.namecol < 0:
        cexit('ERR: microhaps mode needs --namecol option!')

    with open(args.outfile, 'w') as fout:
        mh_name = ''
        mh_seq = ''
        mh_1pos = -1
        mh_2pos = -1
        for entry in positions:
            seq = entry[0]
            pos = int(entry[1])
            name = entry[args.namecol]

            if name == mh_name and mh_seq == seq:
                mh_2pos = pos
                continue

            if mh_name:
                fout.write( '%s\t%d\t%d\t%s\n' % (mh_seq, mh_1pos, mh_2pos, mh_name))

            mh_name = name
            mh_seq = seq
            mh_1pos = pos-1

        fout.write( '%s\t%d\t%d\t%s\n' % (mh_seq, mh_1pos, mh_2pos, mh_name))

    cerr('[I - writing microhap-based BED to %s]' % args.outfile )


