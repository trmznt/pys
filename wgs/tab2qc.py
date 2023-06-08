#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser("Calculate missingness on samples and SNPs")
    p.add_argument('--outprefix',
                   help='Prefix to be used for output files.')
    p.add_argument('infile',
                   help='tabular-format genotype file')

    return p


def tab2qc(args):

    import pandas as pd
    import numpy as np
    from seqpy.core.bioio.tabutils import read_file, write_file, join_metafile

    geno_df = read_file(args.infile)
    alleles = geno_df.geno.get_alleles()
    N, L = alleles.shape

    # count missingness in samples
    sample_qc = np.zeros(N)
    for i in range(N):
        sample_qc[i] = (alleles.iloc[i] == 'X').sum() / L

    sample_qc_df = pd.DataFrame({'Sample': geno_df.geno.get_samples(), 'Missingness': sample_qc})
    sample_qc_df.to_csv(args.outprefix + '-sampleqc.tsv', index=False, sep='\t')

    # count missingness in variants
    variant_qc = np.zeros(L)
    for i in range(L):
        variant_qc[i] = (alleles.iloc[:, i] == 'X').sum() / N

    # prepare variant positions
    pos_df = alleles.columns.to_series().str.split(':', n=1, expand=True)
    pos_df.reset_index(drop=True, inplace=True)

    variant_qc_df = pd.DataFrame({'CHROM': pos_df[0], 'POS': pos_df[1], 'Missingness': variant_qc})
    variant_qc_df.to_csv(args.outprefix + '-variantqc.tsv', index=False, sep='\t')

    import IPython; IPython.embed()



def main(args):
    tab2qc(args)

# EOF
