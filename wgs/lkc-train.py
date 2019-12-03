#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import lkest

import pickle

def init_argparser(p = None):

    if p is None:
        p = arg_parser('classify.py - classify samples based on SNP barcodes')


    p = naltparser.init_argparser(p)

    p.add_argument('--profile', default=None, required=True)
    p.add_argument('--includepos', default=None, required=True)
    p.add_argument('--code', default=None, required=True)
    p.add_argument('--remark', default=None, required=True)
    p.add_argument('--replace', default=False, action='store_true')
    p.add_argument('--update', default=False, action='store_true')

    return p


def train( args ):

    # load profiles if exists

    nalt_parser = naltparser.NAltLineParser(args, datatype='nalt')

    nalt_parser.parse_grouping()
    group_keys = nalt_parser.group_parser.group_keys

    region = nalt_parser.parse_whole()
    samples = nalt_parser.parse_samples()

    poslines = [ line.split() for line in open(args.includepos) ][1:]
    region.filter_poslines(poslines, inplace=False, sort_position=False)
    haplotypes = region.haplotypes()

    cerr('[I - fitting for {}]'.format( args.code ))
    classifier = lkest.SNPLikelihoodEstimator(H0=True)
    classifier.fit(haplotypes, group_keys)
    profile = classifier.get_profile()
    profile.positions = region.P
    profile.code = args.code
    profile.remark = args.remark

    try:
        with open(args.profile, 'rb') as f:
            profiles = pickle.load( f )
    except FileNotFoundError:
        profiles = {}
    if args.code in profiles and not args.replace:
        cexit('ERR: cannot replace ')
    profiles[args.code] = profile.to_dict()
    with open(args.profile, 'wb') as f:
        pickle.dump( profiles, f)

    cerr('[I - profiles saved to {}]'.format( args.profile ))


def main( args ):

    train( args )



