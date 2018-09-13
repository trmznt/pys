#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# calculate population distance matrix D_xy (d_xy - 0.5 * (d_x + d_y))
# 


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser

import allel
import numpy as np

def init_argparser(p=None):

    p = grpparser.init_argparser(p)
    p.add_argument('-o', '--outfile', default='outfile.popdist.txt')
    p.add_argument('infile')

    return p


def main( args ):

    dist2popdist( args )


def dist2popdist( args ):

	# read group assignment

	group_parser = grpparser.GroupParser( args )

	# read distance matrix

	infile = open(args.infile)
	samples = next(infile).split()
	groups = group_parser.assign_groups(samples)
	group_keys = sorted(groups.keys())
	n = len(groups)

	M = np.zeros( (n, n) )

	# calculate intra population
	#for i, g in enumerate(group_keys):
	#	d = c = 0
	#	for x,y in combinations( groups[g], 2):
	#		d += D[x,y]
	#		c += 1
	#	M[i,i] = d/c

	# calculate inter population
	for i, j in combinations_with_replacement(range(n), 2):

		d = c = 0
		for x,y in product(groups[ group_keys[i] ], groups[ group_keys[j] ] ):
			d += D[x,y]
			c += 1
		M[i,j] = d/c

	# perform Dxy calculation

	P = np.zeros( (n,n) )
	for i, j in combinations( range(n), 2 ):

		P[i,j] = P[j,i] = M[i,j] - 0.5(M[i,i] + M[j,j])

	# write distance matrix

	with open(args.outfile + '.popdxy.txt') as outfile:
		# write dxy
		outfile.write( '%s\n' % '\t'.join( group_keys ) )

	with open(args.outfile + '.popdist.txt') as outfile:
		# write distance
		outfile.write( '%s\n' % '\t'.join( group_keys ) )

