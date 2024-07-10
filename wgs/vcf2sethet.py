#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser


def init_argparser():

    p = arg_parser('vcf2sethets - set variants to hets based on certain criteria')

    p.add_argument('--min-alt-count', default=-1, type=int, required=True,
                   help="minimum number of either allele reads to be called hets, "
                   "MalariaGEN value is 2 [-1]")
    p.add_argument('--min-alt-ratio', default=-1, type=float, required=True,
                   help="minimum ratio of either alleles reads to be called hets, "
                   "MalariaGEN value is 0.10 [-1]")
    p.add_argument('--set-het-to-ref', default=False, action='store_true',
                   help='set GT to reference alleles for all het alleles')
    p.add_argument('--set-het-to-alt', default=False, action='store_true',
                   help='set GT to alternate alleles for all het alleles')
    p.add_argument('--set-missing-to-ref', default=False, action='store_true',
                   help='set GT to reference alleles for all missing alleles')
    p.add_argument('--set-missing-to-alt', default=False, action='store_true',
                   help='set GT to alternate alleles for all missing alleles')
    p.add_argument('-o', '--outfile', default='-',
                   help="file output [stdout]")
    p.add_argument('--outlog', default=None,
                   help="output log [None]")
    p.add_argument('infile', nargs='?',
                   help="input VCF file (skipping this arg would default to stdin)")

    return p


def vcf2sethets(args):

    # set variant to hets based on minimum number of minor reads and minimum ratio of
    # minor reads
    # usually depth reassignment is needed after running this script to correct
    # for AC and AF info fields
    # note that this script works for multi-allelic variant as well

    import sys
    import numpy as np
    from cyvcf2 import VCF, Writer

    if not args.infile:
        args.infile = '-'

    min_alt_depth = args.min_alt_count
    min_alt_ratio = args.min_alt_ratio

    vcf = VCF(args.infile)
    w = Writer(args.outfile, vcf)

    w.add_to_header("##spcli_vcf2sethetCmdLine= " + " ".join(sys.argv[1:]))

    logs = []

    for v in vcf:

        AD = v.format('AD')

        if len(v.ALT) == 1:

            minor_depths = AD.min(axis=1)
            minor_args = np.argpartition(AD, -2, axis=1)
            major_alleles = minor_args[:, 1]
            minor_alleles = minor_args[:, 0]

        else:
            # we have multiple alleles, need only to take the all minor depths and
            # check if 1st & 2nd minor depths is the same

            minor_part = np.partition(AD, -2, axis=1)
            minor_args = np.argsort(-AD, axis=1)

            # get the major alleles & temp minor alleles
            major_alleles = minor_args[:, 0]
            minor_alleles = minor_args[:, 1]

            # get cumulative minor depths (eg: all minor depths except the last one)
            minor_depths = minor_part[:, :-1].sum(axis=1)

            # check if 1st minor == 2nd minor, and the value is not zero
            equal_minor_depths = ((minor_part[:, -2] == minor_part[:, -3])
                                  & (minor_depths > min_alt_depth))

            if any(equal_minor_depths):

                # in each sample with equivalent minor depths, reassigned minor allele with
                # minor_args[1] which is guarantee to be smaller allele index than minor_args[2]
                # when the depths of those minor alleles are the same (by np.argsort(-AD))
                minor_alleles[equal_minor_depths] = minor_args[equal_minor_depths, 1]

                # get two of the most common allele index
                allele_indexes = (AD < min_alt_depth).sum(axis=0).argsort()
                #cerr(f'Allele indexes: {allele_indexes}')

                # sanity check if reference is major (allele_indexes[0]) or minor
                # (allele_indexes[1])
                if allele_indexes[0] != 0 and allele_indexes[1] != 0:
                    logs.append(
                        f'WARNING: reference allele is not part of major/minor alleles '
                        f'at variant {repr(v)}'
                    )

        # up to this point, there should be:
        # major_alleles
        # minor_alleles (after adjustment)
        # minor_depths (after adjustment)

        minor_ratios = minor_depths / v.gt_depths

        non_hets = (minor_depths < min_alt_depth) | (minor_ratios < min_alt_ratio)

        # set minor_alleles to major_alleles if minor_depths == 0;
        # run bcftools +setGT -- -t q -n . -i "FORMAT/DP<{mindepth}" to set min depth
        # and recalculate parameters
        null_minor_indexes = (minor_depths == 0)
        minor_alleles[null_minor_indexes] = major_alleles[null_minor_indexes]

        # for all that are non-hets, set both alleles to major allele
        for idx in non_hets.nonzero()[0]:
            allele = major_alleles[idx]
            v.genotypes[idx] = [allele, allele, False]

        # for all that are not non-hets, eg. the hets, set their alleles
        for idx in (~non_hets).nonzero()[0]:
            if args.set_het_to_ref:
                v.genotypes[idx] = [0, 0, False]
                continue
            if args.set_het_to_alt:
                v.genotypes[idx] = [1, 1, False]
                continue
            # sort alleles
            if (maj_allele := major_alleles[idx]) > (min_allele := minor_alleles[idx]):
                alleles = [min_allele, maj_allele, False]
            else:
                alleles = [maj_allele, min_allele, False]
            v.genotypes[idx] = alleles

        for idx in (v.gt_depths == 0).nonzero()[0]:
            if args.set_missing_to_ref:
                v.genotypes[idx] = [0, 0, False]
                continue
            if args.set_missing_to_alt:
                v.genotypes[idx] = [1, 1, False]
                continue
            v.genotypes[idx] = [-1, -1, False]

        # reset genotypes
        v.genotypes = v.genotypes

        w.write_record(v)

    vcf.close()
    w.close()

    cerr(f'Processsed VCF file written to {args.outfile}')

    if args.outlog:
        with open(args.outlog, 'w') as outlog:
            outlog.write('\n'.join(logs))
        cerr(f'Log file written to {args.outlog}')


def main(args):
    vcf2sethets(args)

# EOF
