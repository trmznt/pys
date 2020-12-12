#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import grpparser

import time

try:
    import allel
except:
    cexit('ERR: require properly installed scikit-allel!')

def init_argparser():
    p = arg_parser("Create a ped and map file suitable for isoRelate from VCF and metafile")
    p = grpparser.init_argparser( p )
    p.add_argument('-o', '--outprefix', default='outdata')
    p.add_argument('infile')
    return p

def vcf2ped( args ):
    """ create a ped and map file based on vcf and metafile, suitable for isoRelate """

    # open group file
    group_parser = grpparser.GroupParser( args )

    # open VCF file
    cerr('[I: reading VCF...]')
    start_time = time.monotonic()
    vcfset = allel.read_vcf(args.infile,
                fields = ['samples', 'variants/CHROM', 'variants/POS', 'calldata/GT'])
    cerr('[I: read %s site, %s samples in %d secs]' % (len(vcfset['variants/CHROM']),
         len(vcfset['samples']), time.monotonic() - start_time))

    # assign groups
    samples = vcfset['samples']
    group_parser.assign_groups(samples)
    groups = group_parser.group_keys
    #import IPython; IPython.embed()

    # write to PED
    with open(args.outprefix + '.ped', 'w') as outf:
        for i in range(len(samples)):
            outf.write('%s\t%s\t0\t0\t1\t0\t' % (groups[i], samples[i]))
            alleles = []
            for gt in vcfset['calldata/GT'][:,i]:
                allele_1, allele_2 = gt
                #print(allele_1, allele_2)
                if allele_1 == allele_2:
                    if allele_1 == -1:
                        alleles += [0, 0]
                    elif allele_1 == 0:
                        alleles += [1, 1]
                    elif allele_1 == 1:
                        alleles += [2, 2]
                    else:
                        alleles += [1, 1]
                else:
                    alleles += [1, 2]
            outf.write('\t'.join( str(i) for i in alleles))
            outf.write('\n')
            #import IPython; IPython.embed()

    # write to MAP
    with open(args.outprefix + '.map', 'w') as outf:
        last_pos = 0
        curr_chr = None
        for (chrom, pos) in zip( vcfset['variants/CHROM'], vcfset['variants/POS'] ):
            if curr_chr != chrom:
                curr_chr = chrom
                last_pos = 0
            dist = (pos - last_pos) * 1e-6
            last_pos = pos
            outf.write('%s\t%s:%d\t%8.6f\t%d\n' % (chrom, chrom, pos, dist, pos))



def main( args ):
    vcf2ped( args )
