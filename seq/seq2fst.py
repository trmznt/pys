#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, grpparser, multisequence
from seqpy.core.funcs import profiles

import numpy as np
import itertools
import allel

def init_argparser(p=None):

    if p is None:
        p = arg_parser('seq2pi.py - calculate nucleotide diversity (pi)')

    p = grpparser.init_argparser( p )
    p.add_argument('-o', '--outfile', default='outfile.fst.txt')
    p.add_argument('--sitefile', default='')
    p.add_argument('--nantozero', action='store_true', default=False)
    p.add_argument('infile')

    return p


def main( args ):

    seq2fst( args )


def seq2fst( args ):

    # open and read sequence file
    cerr('[I - reading sequence file %s]' % args.infile)
    seqs = load(args.infile)

    # open and read group/meta file using groupfile/metafile if available
    if args.groupfile or args.metafile:
        cerr('[I - reading group information file]')
        group_parser = grpparser.GroupParser( args )
        group_parser.parse()

        group_seqs = {}

        for seq in seqs:
            try:
                grp = group_parser.group_info[seq.label.decode('ASCII')]
            except KeyError:
                cerr('[W - sample %s is not assign to any group]' % seq.label.decode('ASCII'))
                continue
            if grp in group_seqs:
                group_seqs[grp].append( seq )
            else:
                ms = multisequence()
                ms.append( seq )
                group_seqs[grp] = ms
    else:
        cexit('[ERR - seq2fst.py requires group information!]')

    for grp_seq in group_seqs:
        cerr('[I - group %s has %d sample(s)]' % (grp_seq, len(group_seqs[grp_seq])))

    if args.sitefile:
        # perform FST site-wise
        FST_sites = calc_site_fst( group_seqs, args.nantozero )

        with open(args.sitefile, 'w') as fout:
            for (label, mat) in FST_sites:
                fout.write(label)
                fout.write('\t')
                np.savetxt(fout, mat, fmt='%5.4f', delimiter='\t', newline='\t')
                fout.write('\n')

        cerr('[I - site FST written to %s]' % (args.sitefile))
        return


    FST_mat, groups = calc_fst( group_seqs )

    with open(args.outfile, 'w') as fout:
        fout.write('\t'.join( groups ))
        fout.write('\n')
        np.savetxt(fout, FST_mat, fmt='%5.4f', delimiter='\t')


def calc_fst( mseqs ):

    groups = list(mseqs.keys())
    len_grp = len(groups)
    FST_mat = np.zeros( (len_grp, len_grp) )
    allele_counts = count_allele( mseqs)
    for i,j in itertools.combinations(range(len_grp), 2):

        ac1 = allele_counts[ groups[i] ]
        ac2 = allele_counts[ groups[j] ]

        with np.errstate(divide='ignore', invalid='ignore'):
            num, den = allel.hudson_fst( ac1, ac2 )

            FST_mat[i,j] = np.nanmean( num/den ) #np.sum(num) / np.sum(den)
            FST_mat[j,i] = np.nanstd( num/den )

        cout('%5.4f +- %5.4f : %s <> %s' % (FST_mat[i,j], FST_mat[j,i], groups[i], groups[j]))

    return FST_mat, groups

def calc_site_fst( mseqs, nan_to_zero=False ):

    groups = list(mseqs.keys())
    len_grp = len(groups)
    FST_sites = []
    allele_counts = count_allele( mseqs)
    for i,j in itertools.combinations(range(len_grp), 2):

        ac1 = allele_counts[ groups[i] ]
        ac2 = allele_counts[ groups[j] ]

        with np.errstate(divide='ignore', invalid='ignore'):
            num, den = allel.hudson_fst( ac1, ac2 )

            if nan_to_zero:
                # convert nan to zero
                fst = np.nan_to_num(num/den)
            else:
                fst = num/den
            FST_sites.append( ('%s <> %s' % (groups[i], groups[j]), fst ) )

    return FST_sites


def count_allele(mseqs):

    full_mseqs = multisequence()
    for grp in mseqs:
        full_mseqs.extend( mseqs[grp] )

    na_profiles = profiles.na_profile( full_mseqs, additional_ignore = b'X' )
    consensus_seq = na_profiles.consensus(0.1)

    allele_counts = {}
    for grp in mseqs:
        allele_count = np.zeros( (len(consensus_seq), 2) )
        mseq = mseqs[grp]
        for j in range(len(mseq)):
            seq = mseq[j].seq
            for i in range(len(consensus_seq)):
                if seq[i] == b'N':
                    allele_count[i, 0] += 1
                    allele_count[i, 1] += 1
                elif seq[i] == b'X':
                    continue
                elif seq[i] == consensus_seq[i]:
                    allele_count[i, 0] += 2
                else:
                    allele_count[i, 1] += 2

        allele_counts[grp] = allele_count

    return allele_counts


def to_genotype_array(mseqs):

    # genotype array is array of site vs sample vs [0,0]
    # 0 for the major allele, 1 for the minor allele

    # create matrix profile first
    full_mseqs = multisequence()
    indexes = {}
    idx = 0
    for grp in mseqs:
        full_mseqs.extend( mseqs[grp] )
        indexes[grp] = list(range(idx, idx+len(mseqs[grp])))
        idx += len(mseqs[grp])

    na_profiles = profiles.na_profile( full_mseqs, additional_ignore = b'X' )
    consensus_seq = na_profiles.consensus(0.1)

    genotype_array = np.zeros((len(consensus_seq), len(full_mseqs), 2), dtype=int)
    for i, j in itertools.product(range(len(full_mseqs)), range(len(consensus_seq))):
        seq = full_mseqs[0].seq
        if seq[j] == b'N':
            genotype_array[j,i] == [0, 1]
        elif seq[j] == b'X':
            genotype_array[j, i] == [-1, -1]
        elif seq[j] == consensus_seq[j]:
            genotype_array[j, i] == [0, 0]
        else:
            genotype_array[j, i] == [1, 1]

    return genotype_array, indexes