#!/usr/bin/env spcli

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser

from collections import defaultdict

try:
    from matplotlib import pyplot as plt, colors
except ImportError:
    cexit('ERR: require properly installed matplotlib')


def init_argparser():
    p = arg_parser("Create color annotation for individual based on certain value")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--keep', default=None)
    g.add_argument('--replace', default=None)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--valuefile', default=None)
    g.add_argument('--value', default=None)
    p.add_argument('--column', type=int, default=-1)
    p.add_argument('--normalize', default=False, action='store_true')
    p.add_argument('--reverse', default=False, action='store_true')
    p.add_argument('-o', '--outfile', default='outcolor.txt')

    p.add_argument('infile',
                   help='tab-delimited file')

    return p


def main(args):

    anno2color(args)


def anno2color(args):

    indv_vals = {}

    if args.valuefile and args.column > 0:
        with open(args.valuefile) as valuefile:
            header = next(valuefile)
            cerr('Column value: %s' % header.strip().split('\t')[args.column-1])
            for line in valuefile:
                tokens = line.strip().split('\t')
                indv_vals[tokens[0]] = float(tokens[args.column-1])

    elif args.value:
        indv_vals = defaultdict(lambda: args.value)

    keep_d = {}
    replace_d = {}

    if args.keep:
        with open(args.keep) as keepfile:
            for line in keepfile:
                keep_d[line.strip().split()[0]] = True

    elif args.replace:
        with open(args.replace) as replacefile:
            for line in replacefile:
                replace_d[line.strip().split()[0]] = True

    indv = []
    modify = {}

    with open(args.infile) as infile:
        next(infile)

        if len(keep_d) > 0:

            # keeping individual values
            for line in infile:
                tokens = line.strip().split('\t')
                indv.append(tokens)
                indv_code = tokens[0]
                if indv_code not in keep_d:
                    modify[indv_code] = indv_vals[indv_code]

        elif len(replace_d) > 0:

            # modify individual values
            for line in infile:
                tokens = line.strip().split('\t')
                indv.append(tokens)
                indv_code = tokens[0]
                if indv_code in replace_d:
                    modify[indv_code] = indv_vals[indv_code]

        else:
            cexit('ERR: unable to proceed!')

    if args.reverse:
        values = modify.values()
        min_value = min(values)
        max_value = max(values)
        for k in modify:
            modify[k] = max_value + min_value - modify[k]

    if args.normalize:
        values = modify.values()
        min_value = min(values)
        max_value = max(values)
        min_max = (max_value-min_value)
        for k in modify:
            modify[k] = (modify[k] - min_value) / min_max

    cmap = plt.cm.coolwarm

    with open(args.outfile, 'w') as outfile:
        outfile.write('SAMPLE\tCOLOUR\n')
        for tokens in indv:
            if tokens[0] in modify:
                value = modify[tokens[0]]
                if type(value) == float:
                    value = colors.to_hex(cmap(value))
            else:
                value = tokens[1]
            outfile.write('%s\t%s\n' % (tokens[0], value))

    cerr('Writing to %s' % args.outfile)

# EOF
