""" nalt2nalt.py
    create a new n_alt file based on certain filtering/criteria
"""
from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser

def init_argparser():
    parser = arg_parser('Assess SNP and sample quality from nalt file')
    parser = naltparser.init_argparser(parser)
    parser.add_argument('--remove-indv')
    parser.add_argument('--mac', type=int, default=0)

    return parser


def nalt2nalt(args):
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

    nalt_parser = naltparser.NAltLineParser(args, datatype='nalt')

    # preparing filters
    filters = []
    if args.remove_indv:
        filters.append(
            prepare_remove_indv_filter(args.remove_indv, nalt_parser.parse_samples())
        )
    if args.mac > 0:
        filters.append( prepare_mac_filter(args.mac) )

    whole = nalt_parser.parse_whole(10000)

    for pos, n_alt in whole.parse_positions():
        filtered_n_alt = process_filter(filters, pos, n_alt)
        if filtered_n_alt:
            naltout.write(filtered_n_alt)

def prepare_remove_indv_filter(indv_file, samples):
    indvs = []
    with open(indv_file) as indvfile:
        for indv in indvfile:
            indvs.append( indv.strip() )


def main(args):
    nalt2nalt(args)