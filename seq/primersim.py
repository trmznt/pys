#!/usr/bin/env spcli

import argparse
import os
from seqpy import cerr
from numba import njit
from dataclasses import dataclass


def init_argparser(p=None):
    p = p if p else argparse.ArgumentParser(
        description='primersim - simulate multiplex primers in PCR process'
    )

    p.add_argument('-m', '--minlength', type=int, default=50,
                   help='minimum length of amplicons [50]')
    p.add_argument('-M', '--maxlength', type=int, default=1000,
                   help='maximum length of amplicons [1000]')
    p.add_argument('--maxmismatchprop', type=float, default=0.2,
                   help='maximum proportion of mismatch bases for annealing [0.2]')
    p.add_argument('--reversecomplement', default=False, action='store_true',
                   help='check the reverse-complemented of primers as well')
    p.add_argument('--template',
                   help='reference sequence(s) in FASTA format')
    p.add_argument('-o', '--outfile', default='outprimers.tsv',
                   help='output file name')
    p.add_argument('--outfragment', default='',
                   help='output as fragments in FASTA')
    p.add_argument('infile',
                   help='FASTG file containing primer sequences')

    return p


@njit(nogil=True)
def anneal_primer(template_seq, primer_seq, diff_last_5=2, max_mismatch=0.3):
    """ return a list of (begin_pos, end_pos, direction, score, percentage_diff)
    """

    N_r = len(template_seq)
    N_p = len(primer_seq)
    Max_Diff = int(N_p * max_mismatch)

    temp_results = []
    for i in range(N_r-N_p):
        mismatch = 0
        for j in range(N_p):
            if template_seq[i+j] != primer_seq[j]:
                mismatch += 1
        if mismatch <= Max_Diff:
            temp_results.append((i, i + N_p, mismatch, mismatch / N_p))

    return temp_results


@dataclass
class AnnealPosition(object):
    chrom: str
    begin: int
    end: int
    strand: str
    mismatch: float
    score: float
    primer: object


@dataclass
class Amplicon(object):
    chrom: str
    begin: int
    end: int
    length: int
    anneal_upstream: object
    anneal_downstream: object


def anneal_all_primers(primer_seqs, template_seqs, maxmismatchprop=0.1):

    from seqpy.core.funcs import funcs

    all_primers = []
    for primer_seq in primer_seqs:
        cerr(f'[Priming {primer_seq.label}...]')

        for primer, orientation in [
            (primer_seq.seq, '+'),
            (funcs.reverse_complemented(primer_seq), '-')
        ]:
            for template_seq in template_seqs:
                #cerr(f'[ - Priming against {template_seq.label}]')
                anneal_results = anneal_primer(template_seq.seq, primer,
                                               max_mismatch=maxmismatchprop)
                if any(anneal_results):
                    for begin, end, score, mismatch_prop in anneal_results:
                        all_primers.append(
                            AnnealPosition(
                                chrom=template_seq.label,
                                begin=begin,
                                end=end,
                                strand=orientation,
                                mismatch=mismatch_prop,
                                score=score,
                                primer=primer_seq,
                            )
                        )

    return all_primers


def amplify_all_primers(anneal_list, min_length=100, max_length=1000):

    pos_anneals = [p for p in anneal_list if p.strand == '+']
    neg_anneals = [p for p in anneal_list if p.strand == '-']

    # for each pos anneal, check against neg anneal if they are:
    # - in the same chromosome
    # - can extend up to max_length

    amplicons = []

    for p_a in pos_anneals:
        for q_a in neg_anneals:
            # cerr(f'{p_a[0]} >< {q_a[0]}]')
            if p_a.chrom != q_a.chrom:
                continue
            extension = q_a.end - p_a.begin
            #cerr(f'[ {p_a[0]} - {extension}]')
            if min_length <= extension <= max_length:
                amplicons.append(
                    Amplicon(
                        chrom=p_a.chrom,
                        begin=p_a.begin,
                        end=q_a.end,
                        length=extension,
                        anneal_upstream=p_a,
                        anneal_downstream=q_a
                    )
                )

    return amplicons


def do_simulate(primer_seqs, template_seqs, args):

    annealed_primers = anneal_all_primers(primer_seqs, template_seqs, args.maxmismatchprop)
    return amplify_all_primers(annealed_primers, args.minlength, args.maxlength)


def table_amplicons(amplicons):

    import pandas as pd

    df = pd.DataFrame(
        dict(
            chrom=[], length=[], begin=[], end=[],
            primer_upstream=[], begin_upstream=[], end_upstream=[], seq_upstream=[],
            primer_downstream=[], begin_donwstream=[], end_downstream=[], seq_downstream=[]
        )
    )
    for idx, amp in enumerate(amplicons):
        anneal_up = amp.anneal_upstream
        anneal_dw = amp.anneal_downstream
        df.loc[idx] = [
            amp.chrom, amp.length, amp.begin, amp.end,
            anneal_up.primer.label, anneal_up.begin, anneal_up.end, anneal_up.primer.seq.decode('UTF-8'),
            anneal_dw.primer.label, anneal_dw.begin, anneal_dw.end, anneal_dw.primer.seq.decode('UTF-8'),
        ]

    return df


def fragment_amplicons(amplicons, template_seqs):

    from seqpy.core import bioio

    mseqs = bioio.multisequence()

    for idx, amp in enumerate(amplicons):
        seq = template_seqs.get_by_label(amp.chrom)
        fragment_seq = seq.seq[amp.begin:amp.end]
        fragment_label = (f'{amp.chrom}:{amp.begin}-{amp.end}:'
                          f'{amp.anneal_upstream.primer.label}:'
                          f'{amp.anneal_downstream.primer.label}:'
                          f'{amp.length}')

        mseqs.addseq(fragment_label, fragment_seq)

    return mseqs


def primersim(args):

    from seqpy.core import bioio
    from seqpy.core.funcs import funcs

    # read primer sequences
    primer_seqs = bioio.load(args.infile)

    # read template sequences
    template_seqs = bioio.load(args.template)

    # fix template labels
    for template_seq in template_seqs:
        template_seq.label = template_seq.label.split()[0]

    if args.reversecomplement:
        reverse_complemented_primers = []
        for seq in primer_seqs:
            n_seq = seq.clone()
            n_seq.label = seq.label + '-RC'
            n_seq.seq = funcs.reverse_complemented(seq.seq)
            reverse_complemented_primers.append(n_seq)

        primer_seqs.extend(reverse_complemented_primers)

    amplicons = do_simulate(primer_seqs, template_seqs, args)
    df = table_amplicons(amplicons)

    df.to_csv(args.outfile, sep='\t', index=False)
    cerr(f'[Amplicon report written to {args.outfile}]')

    if args.outfragment:
        mseqs = fragment_amplicons(amplicons, template_seqs)
        bioio.save(mseqs, args.outfragment)
        cerr(f'[Amplicon fragments written to {args.outfragment}]')


def main(args):
    primersim(args)

# EOF
