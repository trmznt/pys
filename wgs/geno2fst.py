#!/usr/bin/env spcli

# Hidayat Trimarsanto <anto@eijkman.go.id>
#
# create an FST file with the following format:
# AVG_FST MAX_FST POP1 POP2 ...


from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core.bioio import tabparser

import allel
import numpy as np

def init_argparser(p=None):

    p = tabparser.init_argparser()
    p.add_argument('-o', '--outfile', default='outfile.fst.txt')

    return p


def main( args ):

    geno2fst( args )


def geno2fst( args ):

    lineparser = tabparser.GenotypeLineParser( args )
    lineparser.set_translator(lineparser.diploid_translator)

    cout('Grouping:')
    groups = lineparser.parse_grouping()
    for k in groups:
        cout(' %12s %3d' % (k, len(groups[k])))

    FST = [] # FST indexed by group_keys
    group_keys = sorted(groups.keys())
    cout(group_keys)

    # output to file
    cout('Writing outfile...')
    outfile = open(args.outfile, 'w')

    outfile.write('CHROM\tPOS\tREGION\tMAX\tMEAN\tMEDIAN\tMAF\t%s\n' % '\t'.join(group_keys) )

    idx = 0
    for (posinfo, genolist) in lineparser.parse():

        idx += 1
        genoarray = allel.GenotypeArray( [genolist]  )

        # calculate MAF
        ac = genoarray.count_alleles()
        num = np.min(ac)
        denom = np.sum(ac)
        if num == denom:
            maf = 0
        else:
            maf = np.min(ac)/np.sum(ac)

        # calculate FST per group against other samples

        fst_sites = []
        for g in group_keys:
            ac_g = genoarray.count_alleles(subpop = groups[g])
            ac_ng = genoarray.count_alleles(subpop = list( lineparser.sample_idx - set(groups[g])))
            num, den = allel.stats.hudson_fst(ac_g, ac_ng)
            fst = num[0]/den[0]
            if not (0.0 <= fst <= 1.0):
                fst = 0
            fst_sites.append( fst )

        if idx % 100 == 0:
            cerr('I: writing position no %d' % idx)

        outfile.write('%s\t%s\t%s\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%s\n' %
                        (posinfo[0], posinfo[1], posinfo[4], np.max(fst_sites), np.mean(fst_sites), np.median(fst_sites), maf,
                            '\t'.join( '%5.4f' % x for x in fst_sites)))

