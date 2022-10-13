#!/usr/bin/env spcli

import argparse
import numpy as np
import pandas as pd

from seqpy import cerr
from seqpy.core.bioio import tabutils

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from matplotlib import pyplot as plt


def init_argparser():

    p = argparse.ArgumentParser()

    p.add_argument('-o', '--outfile')
    p.add_argument('-c', '--colorfile')
    p.add_argument('--outplot')
    p.add_argument('infile')

    return p


def dist2pcoa(args):

    cerr(f'Reading distance matrix from {args.infile}')

    distm_df = tabutils.read_file(args.infile)
    distm_scaled = StandardScaler().fit_transform(distm_df.values)

    pca = PCA(n_components=3)
    pca_feats = pca.fit_transform(distm_scaled)
    pca_df = pd.DataFrame(data=pca_feats, columns=['PC1', 'PC2', 'PC3'])

    if args.colorfile:
        color_df = tabutils.read_file(args.colorfile)
        pca_df['COLOR'] = color_df['COLOR']

    import seaborn as sns

    sns.scatterplot(x='PC1', y='PC2', data=pca_df,
                    c=pca_df['COLOR'])

    plt.savefig(args.outplot)


def main(args):
    dist2pcoa(args)

# EOF
