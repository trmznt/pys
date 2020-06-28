#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser


try:
	from matplotlib import pyplot as plt
except:
	cexit('ERR: require properly installed matplotlib')

import numpy as np

def init_argparser():
	p = arg_parser("Create depth plot")
	p.add_argument('-t', '--title', default='depth')
	p.add_argument('-o', '--outfile', default='outplot.png')

	p.add_argument('infile')

	return p

def main( args ):

	depthplot( args )


def depthplot( args ):

	# read data
	depth_list = []
	with open(args.infile) as fin:
		for line in fin:
			tokens = line.split()
			depth_list.append( (int(tokens[1]), int(tokens[2])) )

	# plot data
	length = depth_list[-1][0] + 1
	x = np.arange(length)
	y = np.zeros(length)

	for (pos, depth) in depth_list:
		y[pos] = depth

	fig, ax = plt.subplots(figsize=(25,5))
	ax.fill_between(x, y, facecolor="lightgreen", color="darkgreen", alpha=0.5)
	ax.set_yscale("log", basey=10)
	ax.set_title(args.title)
	fig.savefig(args.outfile)


