#!/usr/bin/env spcli

import argparse
import numpy as np
import pandas as pd

from seqpy import cerr
from seqpy.core.bioio import tabutils
from numba import njit
from numba_progress import ProgressBar


def init_argparser():

    p = argparse.ArgumentParser()

    p.add_argument('-o', '--outfile')
    p.add_argument('infile')
    return p


@njit(nogil=True)
def njit_loop_calc_dist(allele_m, dist_m, n, allele_length, progress_proxy):
    for i in range(n):
        r_i = allele_m[i]
        for j in range(i - 1):
            d = 0
            c = 0
            r_j = allele_m[j]
            for k in range(allele_length):
                x_k = r_i[k]
                y_k = r_j[k]
                if x_k == 'X' or y_k == 'X':
                    continue
                c += 1
                if x_k != y_k:
                    d += 1
            dist_m[i, j] = dist_m[j, i] = (d / c)
            progress_proxy.update(1)


def calc_dist(allele_df):
    allele_m = allele_df.values.astype(str)
    n = len(allele_m)
    dist_m = np.zeros((n, n))
    variant_length = len(allele_m[0])
    total_iter = n * (n - 1) / 2
    cerr(f'Calculating distance for {n} samples with {variant_length} variants...')

    with ProgressBar(total=total_iter) as pbar:
        njit_loop_calc_dist(allele_m, dist_m, n, variant_length, pbar)

    return dist_m


def tab2dist(args):

    cerr(f'Reading genomic data from {args.infile}')
    df = tabutils.read_file(args.infile)
    allele_df = df.geno.get_alleles()
    sample_df = df.geno.get_samples()

    distm = calc_dist(allele_df)

    distm_df = pd.DataFrame(data=distm, columns=sample_df)
    tabutils.write_file(args.outfile, distm_df)
    cerr(f'Output is written to {args.outfile}')


def main(args):
    tab2dist(args)

# EOF
