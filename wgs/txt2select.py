#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser

try:
    import pandas
except:
    cexit('ERR: require properly installed pandas')
import numpy as np

def init_argparser():
    p = arg_parser("Select CHROM and POS based on some criteria")
    p.add_argument('--max', type=int, default=100)
    p.add_argument('--column', required=True)
    p.add_argument('-o', '--outfile', default="outfile.pos.txt")

    p.add_argument('infile')

    return p


def main( args ):

    txt2select( args )


def txt2select( args ):

	df = pandas.read_table(args.infile, delimiter='\t', na_values='  nan')

	filtered_positions = {}

	for i in args.column.split(','):

		i = int(i)
		df_filtered = df.nlargest(args.max, df.columns[i])
		for r in df_filtered.itertuples():
			filtered_positions[(r[1], r[2])] = True

	with open(args.outfile, 'w') as outfile:
		outfile.write('CHROM\tPOS\n')
		for k in filtered_positions:
			outfile.write('%s\t%s\n' % (k))




