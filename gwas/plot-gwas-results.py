#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('run-fastlmm - run FastLMM')
    p.add_argument('--outqq', default='qq.png')
    p.add_argument('--outmht', default='mht.png')
    p.add_argument('--add-title-line', action='append', default=[])
    p.add_argument('infile')
    return p


def mhtplot(results_df, outplot, title=None):

    import fastlmm.util.util as flutil
    import matplotlib.pyplot as plt

    plt.rcParams['figure.figsize'] = (10.0, 8.0)
    bonferroni_pvalue = 0.05 / len(results_df)
    flutil.manhattan_plot(results_df[["Chr", "ChrPos", "PValue"]].values,
                          pvalue_line=bonferroni_pvalue,
                          xaxis_unit_bp=False)
    if title:
        plt.title(title, loc='left')
    plt.savefig(outplot)
    plt.close()



def estimate_lambda(p_values):
    import numpy as np
    import scipy.stats as st

    LOD2 = np.median(st.chi2.isf(p_values, 1))
    return LOD2 / 0.456


def qqplot(p_values, outplot=None, title=None):
    import numpy as np
    import seaborn as sns
    from matplotlib import pyplot as plt


    p_values = p_values.copy()
    M = len(p_values)
    pnull = (0.5 + np.arange(M)) / M
    p_values[p_values > 1] = 1

    qnull = -np.log10(pnull)
    qemp = -np.log10(np.sort(p_values))

    #ax = sns.scatterplot(x=qnull, y=qemp, s=2)
    g = sns.JointGrid(x=qnull, y=qemp, marginal_ticks=True)
    g.plot_joint(sns.scatterplot, s=2, linewidth=0)
    g.plot_marginals(sns.histplot)
    ax = g.ax_joint
    max_value = max(max(qnull), max(qemp))
    ax.plot([0, max_value], [0, max_value], linewidth=0.5, c='#cecede')
    est_lambda = estimate_lambda(p_values)
    ax.text(0, max(qemp), f'est lambda = {est_lambda:7.6f}')
    ax.set_xlabel('expected -log10(P)')
    ax.set_ylabel('observed -log10(P)')

    plt.suptitle(title, x=0.1, ha='left', fontsize='medium')
    plt.tight_layout()

    if outplot:
        plt.savefig(outplot)
    else:
        plt.show()


def plot_gwas_results(args):

    import pandas as pd

    df = pd.read_feather(args.infile)

    # add title line
    title = '\n'.join(args.add_title_line)

    mhtplot(df, args.outmht, title=title)
    qqplot(df.PValue, args.outqq, title=title)


def main(args):
    plot_gwas_results(args)


# EOF
