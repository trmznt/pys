#!/usr/bin/sh spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser

import numpy as np

import io, time

def init_argparser(p=None):

    p = naltparser.init_argparser(p)

    p.add_argument('-o', '--outfile')
    p.add_argument('--type', default='nalt', choices=['nalt', 'ralt'])
    return p


def ralt2pickle( args ):

    alt_parser = naltparser.NAltLineParser( args, datatype=args.type
        , with_group=False, with_position=False )

    whole_region = alt_parser.parse_whole()
    whole_region.df_M.to_pickle(args.outfile)

def main(args):

    ralt2pickle( args )
