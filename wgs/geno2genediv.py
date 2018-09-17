#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an FST file with the following format:
# CHROM POS GENE SNPS HAPLOTYPES FST AVG_FST MAX_FST MIN_FST POP1 POP2 ...


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.genediv.txt')

    return p


def main( args ):

    geno2genediv( args )


def geno2genediv( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.diploid_translator)

    # set group
    groups = lineparser.parse_grouping()

    cout('Grouping:')
    group_keys = sorted(groups.keys())
    for k in group_keys:
        cout(' %12s %3d' % (k, len(groups[k])))

    outfile = open(args.outfile, 'wt')
    outfile.write('CHROM\tPOS\tREGION\tN_SNP\tN_HAPLO\tMEAN\tMEDIAN\tMAX\tMIN\t%s\n' %
    				'\t'.join( group_keys ))

    for idx, region in enumerate(lineparser.parse_genes()):
    	haplotypes = set( region.haplotypes())
    	enc_haplos = region.encode_haplotypes()
    	assert len(haplotypes) == max(enc_haplos) + 1
    	haploarray = allel.HaplotypeArray( [enc_haplos] )

    	cerr( 'I: calculating %d - %s' % (idx, region.name))

    	value = []
    	for g in group_keys:
            ac_g = haploarray.count_alleles(subpop = groups[g])
            ac_ng = haploarray.count_alleles(subpop = list( lineparser.sample_idx - set(groups[g])))
            num, den = allel.stats.hudson_fst(ac_g, ac_ng)
            value.append(den)

    	#print(idx, '%4d' % len(haplotypes), max(enc_haplos), region.name, value)
    	params = ( np.mean(value), np.median(value), np.max(value), np.min(value))
    	outfile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\n' % (
    			region.P[0][0], region.P[0][1], region.name, len(region.P), len(haplotypes),
    			'\t'.join( '%5.4f' % x for x in params),
    			'\t'.join( '%5.4f' % x for x in value)))

