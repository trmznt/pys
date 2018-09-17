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
    p.add_argument('-o', '--outfile', default='outfile.genehe.txt')

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
    outfile.write('CHROM\tPOS\tREGION\tN_SNP\tN_HAPLO\tFST\tdHe\tHe\tMEAN\tMEDIAN\tMAX\tMIN\t%s\n' %
    				'\t'.join( group_keys ))

    for idx, region in enumerate(lineparser.parse_genes()):
        haplotypes = set( region.haplotypes())
        enc_haplos = region.encode_haplotypes()
        haploarray = allel.HaplotypeArray( [enc_haplos] )

        cerr( 'I: calculating %d - %s' % (idx, region.name))

        # calculate total He first
        He = 1 - np.sum( haploarray.count_alleles().to_frequencies()**2 )

        # calculate He per population, He_p
        values = []
        pHe = 0
        for g in group_keys:

            he_p = 1 - np.sum(
                haploarray.count_alleles(subpop=groups[g]).to_frequencies()**2 )
            pHe += he_p * len(groups[g])
            values.append(he_p)

        dHe = He - pHe / sum( len(x) for x in groups.values() )
        FST = dHe/He

    	#print(idx, '%4d' % len(haplotypes), max(enc_haplos), region.name, value)
        params = ( FST, dHe, He, np.mean(values), np.median(values), np.max(values), np.min(values))
        outfile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\n' % (
                region.P[0][0], region.P[0][1], region.name, len(region.P), len(haplotypes),
                '\t'.join( '%5.4f' % x for x in params),
                '\t'.join( '%5.4f' % x for x in values)))
