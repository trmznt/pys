#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser


def init_argparser():
    p = arg_parser("Create colour annotation based on header of a file")
    p = grpparser.init_argparser(p)
    p.add_argument('--delimiter', default='\t')
    p.add_argument('-o', '--outfile', default="outfile.anno")

    p.add_argument('infile')

    return p


def main( args ):

    grp2anno( args )


def grp2anno( args ):

    # read group file

    group_parser = grpparser.GroupParser( args )

    # open infile
    with open(args.infile) as infile:
        header = next(infile)

    if args.delimiter is not None:
        samples = header.strip().split(args.delimiter)
    else:
        samples = header.strip().split()

    #import IPython; IPython.embed()

    groups = group_parser.assign_groups(samples)
    group_keys = sorted(groups.keys())
    colours = group_parser.colour_list()

    with open(args.outfile + '.indv.txt', 'w') as outfile:
        outfile.write('SAMPLE\tCOLOUR\n')
        for s,c in zip(samples, colours):
            outfile.write('%s\t%s\n' % (s,c))

    with open(args.outfile + '.group.txt', 'w') as outfile:
        outfile.write('GROUP\tCOLOUR\n')
        for g, c in zip( group_keys, group_parser.group_colour_list(group_keys)):
            outfile.write('%s\t%s\n' % (g, c))
