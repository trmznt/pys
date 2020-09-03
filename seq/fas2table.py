#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, grpparser, multisequence


def init_argparser(p=None):

    if p is None:
        p = arg_parser('fas2table.py - generate mutation table')

    p.add_argument('-o', '--outfile', default='outfile.txt')
    p.add_argument('--reffile', required=True)
    p.add_argument('infile')

    return p


def main( args ):

    fas2table( args )


def fas2table(args):

	msa = load(args.infile)
	ref = load(args.reffile)

	table = generate_table(msa, ref)

	with open(args.outfile, 'w') as fout:
		for (label, muts) in table:
			fout.write('%s/\t%s\n' % (label, ' '.join(muts)))

	cerr('[Writing table to %s]' % args.outfile) 


def generate_table(msa, ref):

	refseq = ref[0].seq.upper()
	X = ord('X')

	table = []

	for seq in msa:

		s = seq.seq.upper()
		muts = []
		for i, p in enumerate(zip(refseq, s), 1):
			if p[1] == X: continue
			if p[0] != p[1]:
				muts.append( '%s%d%s' % (chr(p[0]), i, chr(p[1])))
		table.append( (seq.label, muts) )

	return table


