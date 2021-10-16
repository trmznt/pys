#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outplot', default='outplot.pdf')
    p.add_argument('--mindepth', type=int, default=0)
    p.add_argument('--debug', default=False, action='store_true')
    p.add_argument('infile')
    return p


def generate_heatmap(data, title, outfile):

    # calculate the size for the plot
    L = len(data.columns)
    N = len(data.index)
    max_ylabel_len = data.index.str.len().max()
    max_xlabel_len = data.columns.str.len().max()
    plt.subplots(figsize=((L * 2 + max_ylabel_len / 2) / 10, (N * 2 + max_xlabel_len / 2) / 10))
    ax = sns.heatmap(data, yticklabels=1, xticklabels=1, square=True, norm=LogNorm(),
                     cmap='Wistia', linewidths=0.25, cbar_kws={'shrink': 0.5})
    # ax.set_xticklabels(data.columns, fontsize=2)
    # ax.set_yticklabels(data.index, fontsize=2)
    ax.set_title(title)
    plt.tight_layout()
    # plt.show()
    plt.savefig(outfile, dpi=300)
    plt.close()


def depths2heatmap(args):

    snp_depths = pd.read_table(args.infile, sep='\t', index_col=0)
    title = 'Depth Plot of SNPs and Regions'
    if args.mindepth > 0:
        snp_depths.values[snp_depths.values < args.mindepth] = 0
        title = title + f' (> {args.mindepth})'
    generate_heatmap(snp_depths, title, args.outplot)


def main(args):
    depths2heatmap(args)


# this script can be run independently or through seqpy spcli
if __name__ == '__main__':
    args = init_argparser().parse_args()
    if args.debug:
        from ipdb import launch_ipdb_on_exception
        with launch_ipdb_on_exception:
            main(args)
    else:
        main(args)

# EOF
