#!/usr/bin/env spcli

import argparse
from seqpy import cerr


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('--max-alt-alleles', type=int, default=8)
    p.add_argument('--sample-column', default=False, action='store_true')
    p.add_argument('-o', '--outfile')
    p.add_argument('infile')

    return p


def zarr2depthtable(args):

    import pandas as pd
    from seqpy.core.sgk import sgio

    ds = sgio.load_dataset(args.infile, args.max_alt_alleles)

    if not args.sample_column:
        df = pd.DataFrame(
            ds.call_DP.values.T,
            columns=[
                f'{c}:{p}' for c, p in zip(ds.contig_id.values[ds.variant_contig],
                                           ds.variant_position.values)
            ]
        )

        # add sample column
        df.insert(0, 'SAMPLE', ds.sample_id.values)

    else:
        df = pd.DataFrame(
            ds.call_DP.values,
            columns=ds.sample_id.values
        )

        # add variant/marker solumn
        df.insert(
            0,
            'MARKER',
            [
                f'{c}:{p}' for c, p in zip(ds.contig_id.values[ds.variant_contig],
                                           ds.variant_position.values)
            ])

    df.to_csv(args.outfile, index=False, sep='\t')
    cerr(f'[Depth table written to {args.outfile}]')


def main(args):
    zarr2depthtable(args)


# EOF
