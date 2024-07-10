#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('prepare-phenotype - prepare phenotype & covariate variables')
    p.add_argument('--samplefile')
    p.add_argument('--phenotype')
    p.add_argument('--covariate', default=[], action='append')

    # phenotype transformation
    p.add_argument('--to-ordinal', type=int, default=-1)

    # output  options
    p.add_argument('--outsample')
    p.add_argument('--outcovars')
    p.add_argument('--outphenotype')

    p.add_argument('infile')
    return p


def prepare_phenotype(args):

    import pandas as pd
    import numpy as np

    # prepare covariate + phenotype
    covars_df = pd.read_feather(args.infile)
    covars_df = covars_df.loc[:, ['WGSID'] + [args.phenotype] + args.covariate]

    sample_df = pd.read_table(args.sample, header=None)
    sample_df.rename(columns={0: 'WGSID'}, inplace=True)
    merged_df = sample_df.merge(covars_df)
    #import IPython; IPython.embed()
    #merged_df['log_CQ_EC50'] = np.log(merged_df.CQ_EC50)
    #merged_df.loc[merged_df.CQ_EC50 > 300, 'CQ_EC50'] = 300
    merged_df = merged_df.dropna()

    # write phenotype
    if args.outphenotype:

        if args.to_ordinal > 0:
            # change phenotype value to ordinal using KMeans
            phenotype_df['transformedPhenotype'] = None
            args.phenotype

        phenotype_df = merged_df.loc[:, ['WGSID', 'WGSID'] + [args.phenotype]]
        phenotype_df.to_csv(args.outphenotype, header=False, index=False, sep='\t')

    # write covariates
    if args.outcovars:
        covariates_df = merged_df.loc[:, ['WGSID', 'WGSID'] + args.covariate]
        covariates_df.to_csv(args.outcovars, header=False, index=False, sep='\t')

    # write samples
    #with open(output.sample_fn, 'wt') as sample_out:
    #    sample_out.write('#IID\n')
    #    sample_out.write('\n'.join(merged_df.WGSID))
    if args.outsample:
        merged_df['#FID'] = merged_df.WGSID
        merged_df['IID'] = merged_df.WGSID
        merged_df[['#FID', 'IID']].to_csv(args.outsample, index=False, sep='\t')


def main(args):
    prepare_phenotype(args)


# EOF
