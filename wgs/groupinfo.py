#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser

def init_argparser(p = None):

    if p is None:
        p = arg_parser('groupinfo.py - provide group information')

    p = grpparser.init_argparser( p )
    p.add_argument('--fmt', default='text', choices=['pickle', 'text', 'list'])
    p.add_argument('infile')
    return p


def groupinfo( args ):

    # open and read the first line of infile
    if args.fmt == 'pickle':

        from seqpy.core.bioio import naltparser
        from types import SimpleNamespace
        nalt_args = SimpleNamespace( infile = args.infile, fmt = 'pickle', n = -1)
        nalt_parser = naltparser.NAltLineParser(nalt_args
                , with_group=False, with_position=False)
        samples = nalt_parser.samples

    elif args.fmt == 'list':
        with gzopen(args.infile) as f:
            buf = f.read()
            samples = buf.split()

    else:
        with gzopen(args.infile) as f:
            samples = f.readline().strip().split()

    group_parser = grpparser.GroupParser( args )
    groups = group_parser.assign_groups(samples)
    total = 0
    cout('Groups:')
    for g in sorted(groups.keys()):
        c = len(groups[g])
        cout('  %3d - %s' % (c, g))
        total += c
    cout('Total: %d samples' % total)

def main( args ):
    groupinfo( args )