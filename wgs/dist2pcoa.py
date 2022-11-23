
from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument('-o', '--outfile')
    p.add_argument('-n', '--n_components', type=int, default=3)
    p.add_argument('-m', '--metafile')
    p.add_argument('infile')
    return p


def dist2pcoa(args):

    from seqpy.core.bioio.tabutils import read_file, write_file, join_metafile
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import pandas as pd

    cerr(f'[Reading distance matrix from {args.infile}')

    distm_df = read_file(args.infile)
    samples = distm_df.columns

    # prepare metadata
    if args.metafile:
        sample_df, errs = join_metafile(samples, args.metafile)

    # calculate PCA
    scaled_distm = StandardScaler().fit_transform(distm_df.values)
    pca = PCA(n_components=args.n_components)
    pca_feats = pca.fit_transform(scaled_distm)
    pca_df = pd.DataFrame(
        data=pca_feats,
        columns=[f'PC{x + 1}' for x in range(args.n_components)]
    )

    # join sample_df and pca_df
    df = pd.concat([sample_df, pca_df], axis=1)
    write_file(args.outfile, df)
    cerr(f'[PCoA written to {args.outfile}]')


def main(args):
    dist2pcoa(args)

# EOF
