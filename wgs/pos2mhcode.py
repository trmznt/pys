
from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile')
    p.add_argument('infile')
    return p


def pos2mhcode(args):

    from seqpy.core.bioio.posutils import read_posfile
    from seqpy.core.bioio.tabutils import write_file

    pos_df = read_posfile(args.infile)
    mhcode_df = pos_df.pos.to_mhcode()
    write_file(args.outfile, mhcode_df)
    cerr(f'[microhaplotype codes written to {args.outfile}]')


def main(args):
    pos2mhcode(args)

# EOF
