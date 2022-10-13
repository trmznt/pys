#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser
from seqpy.core.bioio import metautils, tabutils


def init_argparser():
    p = arg_parser("Create colour annotation based on header of a file")
    p.add_argument('--metafile', required=True,
                   help='TSV/CSV file containing SAMPLES and group columns')
    p.add_argument('--specfile')
    p.add_argument('-s', default=False, action='store_true',
                   help='include sample size of each group in legend')
    p.add_argument('-o', '--outfile', default="outfile.anno")

    p.add_argument('infile')

    return p


def tab2anno(args):

    # read infile using tabutils
    df = tabutils.read_file(args.infile)

    # use geno extension, if error then just read headers
    try:
        samples = df.geno.get_samples()
    except AttributeError:
        samples = df.columns
    cerr(f'Reading {len(samples)} samples from {args.infile}')

    # if read the metadata file
    # get filename and specs
    metafile, columns = args.metafile.split(':')
    meta_df = tabutils.read_file(metafile)
    if not columns:
        columns = meta_df.columns[[0, 1]]
    else:
        columns = columns.split(',')

    joined_df = meta_df.meta.join_to_samples(samples, columns)
    cerr(f'Obtaining {len(joined_df)} samples after joined with metadata '
         f'{metafile}')

    # get spec file and its specs
    specfile, column = args.specfile.split(':')
    if not column:
        raise ValueError(f'Please provide the column key at the end of filename, eg: '
                         f'{specfile}:COLUMN_NAME')
    spec_df = tabutils.read_file(specfile)
    joined_df = joined_df.meta.join(spec_df, column)
    cerr(f'Obtaining {len(joined_df)} samples after joined with spec file '
         f'{specfile}')

    tabutils.write_file(args.outfile, joined_df)
    cerr(f'Output written to {args.outfile}')


def main(args):
    tab2anno(args)

# EOF
