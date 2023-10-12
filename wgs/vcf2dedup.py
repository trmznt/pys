#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser('vcf2dedup - removing duplicated positions from a VCF file')

    p.add_argument('--keep-max-mac', default=False, action='store_true',
                   help="keep the position with highest MAC")
    p.add_argument('-o', '--outfile', default='-',
                   help="file output [stdout]")
    p.add_argument('infile',
                   help="input VCF file, use - for piped stdin")

    return p


class SimpleCounter(object):

    def __init__(self, start=0):
        self.counter = start

    def add(self, value=1):
        self.counter += value

    def __str__(self):
        return str(self.counter)


def vcf2dedup(args):

    import sys
    from cyvcf2 import VCF, Writer

    vcf = VCF(args.infile)
    w = Writer(args.outfile, vcf)

    w.add_to_header("##spcli_vcf2dedupCmdLine= " + " ".join(sys.argv[1:]))

    dedup_count = SimpleCounter()

    if args.keep_max_mac:

        curr_chrom = None
        prev_position = -1
        duplicated_variants = []

        def _check_duplicated_variant():
            if len(duplicated_variants) == 1:
                w.write_record(duplicated_variants[0])
                duplicated_variants.clear()

            else:
                kept_variant = duplicated_variants[0]
                for prev_variant in duplicated_variants[1:]:
                    if kept_variant.INFO.get('AC') < prev_variant.INFO.get('AC'):
                        kept_variant = prev_variant
                w.write_record(kept_variant)
                duplicated_variants.clear()
                dedup_count.add()

        for v in vcf:

            if curr_chrom != v.CHROM:
                curr_chrom = v.CHROM
                prev_position = -1

            if prev_position > v.start:
                cexit('[Fatal Error - VCF file is not sorted by position]')

            if prev_position == v.start:
                duplicated_variants.append(v)
                continue

            # after this line, we processes the previous line
            if any(duplicated_variants):
                _check_duplicated_variant()
            duplicated_variants.append(v)
            prev_position = v.start

        if any(duplicated_variants):
            _check_duplicated_variant()

    else:
        cexit('[Fatal Error - please indicate which position to be retain for'
              'duplicated positions]')

    cerr(f'[Deduplicated {dedup_count} positions]')
    vcf.close()
    w.close()


def main(args):
    vcf2dedup(args)

# EOF



