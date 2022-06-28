#!/usr/bin/env spcli

import argparse
import os

from seqpy import cerr
from seqpy.core.sgk import sgio, sgutils
from seqpy.core.bioio import posutils, tabutils


def init_argparser():
    p = argparse.ArgumentParser()
    p.add_argument('infile')

    return p


def zarr2info(args):

    ds = sgio.load_dataset(args.infile)
    for k in ds.dims:
        cerr(f'    {k}: {ds.dims[k]}')


def main(args):
    zarr2info(args)

# EOF
