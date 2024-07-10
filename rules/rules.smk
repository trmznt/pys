

# get params from config

known_variant_dir = config.get('known_variant_dir')
refseq = config.get('refseq', "")

REGIONS = config.get('regions', [])
complete_region = config.get('complete_region', 'complete')
# use empty list to indicate unused variables
variant_file = config.get('variant_file', [])
sample_file = config.get('sample_file', [])
strict_samples = config.get('strict_samples', True)
meta_file = config.get('metadata_file', [])
chromtranslation_file = config.get('chromtranslation_file', 'NOFILE')
target_population = config.get('target_population', '')

#MAC = config.get('MAC', 3)
#source_dir = config.get('source_dir', None)
#sample_threshold = float(config.get('sample_threshold', 0.50))
#variant_threshold = float(config.get('variant_threshold', 0.50))
mindepth = int(config.get('mindepth', 5))

alphanum = r'\w+'
alphanumdash = r'[\w-]+'
alphanumdashdotplus = r'[\.+\w-]+'
digit = r'\d+'  # integers
digitdot = r'[\.\d]+'   # float, fraction

wildcard_constraints:
    reg = alphanumdash,        # chromosome region name
    complete_region = alphanumdash,
    population=alphanumdashdotplus,
    fname = r'[\.+\w-]+',   # standard file name
    MAC = digit,
    MAF = digitdot,
    mindepth = r'\d+',
    st = r'[\.\d]+',
    vt = r'[\.\d]+',
    qual = digitdot,
    fws = digitdot,
    var = r'[\w]+',         # variable name
    grp = alphanum,
    dhets = digit,
    rhets = digitdot,


color_palettes = {
    "xgfs_normal12": [
        "#ebac23",
        "#b80058",
        "#008cf9",
        "#006e00",
        "#00bbad",
        "#d163e6",
        "#b24502",
        "#ff9287",
        "#5954d6",
        "#00c6f8",
        "#878500",
        "#00a76c",
        "#bdbdbd"
      ]
}
# for color palette, see https://gist.github.com/xgfs/37436865b6616eebd09146007fea6c09


# note also some directives and built-in functions:
# workflow.basedir
# workflow.snakefile
# srcdir()

# -- sample and variant handling --


rule flt_pass_depth:
    threads: 3
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        idx = "{pfx}/vcfs/{reg}.vcf.gz.tbi"
    output:
        '{pfx}-PASS-d{mindepth}/vcfs/{reg}.vcf.gz'
    shell:
        'bcftools view -R {known_variant_dir}/{wildcards.reg}.bed.gz {input.vcf} '
        '| bcftools +setGT -- -t q -n . -i "FORMAT/DP<{wildcards.mindepth}" '
        '| bcftools view -o {output}'


rule flt_pass:
    threads: 2
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        idx = "{pfx}/vcfs/{reg}.vcf.gz.tbi"
    output:
        '{pfx}-PASS/vcfs/{reg}.vcf.gz'
    shell:
        'bcftools view -R {known_variant_dir}/{wildcards.reg}.bed.gz {input.vcf} '
        '| bcftools view -o {output}'


rule set_hets_depth:
    threads: 2
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        idx = "{pfx}/vcfs/{reg}.vcf.gz.tbi"
    output:
        '{pfx}-HETd{dhets}r{rhets}-d{mindepth}/vcfs/{reg}.vcf.gz'
    shell:
        'spcli $PYS/wgs/vcf2sethet.py --min-alt-count {wildcards.dhets} '
        '--min-alt-ratio {wildcards.rhets} {input.vcf}'
        '| bcftools +setGT -- -t q -n . -i "FORMAT/DP<{wildcards.mindepth}" '
        '| bcftools view -o {output}'


rule set_depth:
    threads: 2
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        idx = "{pfx}/vcfs/{reg}.vcf.gz.tbi"
    output:
        '{pfx}-d{mindepth}/vcfs/{reg}.vcf.gz'
    shell:
        'bcftools +setGT {input.vcf} -- -t q -n . -i "FORMAT/DP<{wildcards.mindepth}" '
        '| bcftools view -o {output}'


ruleorder: flt_pass_depth > flt_pass > set_hets_depth > set_depth


rule flt_varqual:
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        '{pfx}-Q{qual}/vcfs/{reg}.vcf.gz'
    shell:
        "bcftools view -e 'QUAL < {wildcards.qual}' -o {output} {input}"


rule flt_MAC:
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-MAC{MAC}/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -e "MAC<{wildcards.MAC}" -o {output} {input}'


rule flt_MAF:
    # this rule filters variants based on minimum Minor Allele Frequency (MAF)
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-MAF{MAF}/vcfs/{reg}.vcf.gz"
    shell:
        # 'bcftools +fill-tags {input} -- -t AF '
        # '| bcftools view -e "MAF<{wildcards.MAF}" -o {output} {input}'
        'bcftools +fill-tags {input} -- -t AF '
        '| bcftools view -q {wildcards.MAF}:minor  -o {output}'


rule flt_biallelic:
    # this rule filters for bi-allelic variants
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-biallelic/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -m2 -M2 -o {output} {input}'


rule vcf_atomize:
    # this rule splits complex variants and multi-allele SNVs to individual
    # VCF lines, adding TYPE and F_MISSING to INFO field
    threads: 2
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        ref = refseq
    output:
        "{pfx}-atom/vcfs/{reg}.vcf.gz"
    log:
        log1 = "{pfx}-atom/logs/norm-a-m-snps-f-{reg}.log",
        log2 = "{pfx}-atom/logs/fill-tags-TYPE-F_MISSING-{reg}.log"
    shell:
        "bcftools norm -a -m -snps -f {input.ref} {input.vcf} 2> {log.log1} "
        "| bcftools +fill-tags -o {output} -- -t TYPE,F_MISSING 2> {log.log2}"


rule vcf_dedup:
    # this rule remove duplicated positions and keeping the variant that
    # has bigger minor allele count
    # combination of -atom-dedup is identical to mark indels and other alternate
    # alleles (eg. non-biallelic) as missing values
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-dedup/vcfs/{reg}.vcf.gz"
    log:
        "{pfx}-dedup/logs/dedup-{reg}.log"
    shell:
        "gunzip -c {input} | spcli $PYS/wgs/vcf2dedup.py --keep-max-mac - 2> {log} "
        "| bgzip -c > {output}"


rule flt_noindels:
    # this rule removes indels and other non-SNV variants
    threads: 2
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-noindels/vcfs/{reg}.vcf.gz"
    log:
        "{pfx}-noindels/logs/view-V-indels-{reg}.log"
    shell:
        "bcftools view -V indels,other -o {output} {input} 2> {log}"


rule flt_snv:
    # this rule filter-in only single nucleotide variants
    threads: 2
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-snv/vcfs/{reg}.vcf.gz"
    log:
        "{pfx}-snv/logs/view-v-snps-{reg}.log"
    shell:
        "bcftools view -v snps -o {output} {input} 2> {log}"


rule flt_normalize_snp:
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        "{pfx}-normsnv/vcfs/{reg}.vcf.gz"
    shell:
        "bcftools norm -a {input} "
        "| bcftools view -v snps -o {output}"


rule flt_variants_samples:
    threads: 3
    input:
        vcf = "{src_dir}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{src_dir}/vcfs/{reg}.vcf.gz.csi",
        variant_file = variant_file,
        sample_file = sample_file,
    output:
        "{src_dir}-variants-samples/vcfs/{reg}.vcf.gz"
    log:
        log1 = "{src_dir}-variants-samples/logs/bcftools-view-R-{reg}.log",
        log2 = "{src_dir}-variants-samples/logs/bcftools-view-S-{reg}.log"
    shell:
        'bcftools view -R {input.variant_file} {input.vcf} 2> {log.log1} '
        '| bcftools view -S {input.sample_file} -o {output} 2> {log.log2}'


rule flt_samples:
    # filter VCFs by samples
    # requires:
    #   source_dir
    #   samples_file [a text file containing sample ids]
    threads: 2
    input:
        vcf = "{src_dir}/vcfs/{reg}.vcf.gz",
        sample_file = sample_file
    output:
        "{src_dir}-samples/vcfs/{reg}.vcf.gz"
    params:
        strict_flag = '' if strict_samples else '--force-samples'
    shell:
        'bcftools view -S {input.sample_file} {params.strict_flag} -o {output} {input.vcf}'


rule flt_variants:
    # filter VCFs by variants
    threads: 2
    input:
        vcf = "{src_dir}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{src_dir}/vcfs/{reg}.vcf.gz.tbi",
        variant_file = variant_file
    output:
        "{src_dir}-variants/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -R {input.variant_file} -o {output} {input.vcf}'


ruleorder: flt_variants_samples > flt_variants > flt_samples


rule flt_V:
    # filter variants based on missingness threshold
    threads: 3
    input:
        variant = "{pfx}/variants_{vt}.txt",
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        "{pfx}-V{vt}/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -R {input.variant} {input.vcf} '
        '| bcftools view -e "F_MISSING>{wildcards.vt}" -o {output} '


rule prep_V:
    # prepare variant list for certain threshold
    localrule: True
    input:
        variant = "{pfx}/QC/variants.tsv"
    output:
        variant = "{pfx}/variants_{vt}.txt"
    run:
        import pandas as pd
        df_variant = pd.read_table(input.variant)
        df_variant[df_variant.F_MISS <= float(wildcards.vt)][['CHROM', 'POS']].to_csv(
            output.variant, sep='\t', index=False, header=False
        )


rule flt_S:
    # filter samples based on missingness threshold
    threads: 3
    input:
        sample = "{pfx}/samples_{st}.txt",
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        "{pfx}-S{st}/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -S {input.sample} -o {output} {input.vcf}'


rule prep_S:
    # prepare sample list for certain threshold
    localrule: True
    input:
        sample = "{pfx}/QC/samples.tsv"
    output:
        sample = "{pfx}/samples_{st}.txt"
    run:
        import pandas as pd
        df_sample = pd.read_table(input.sample)
        df_sample[df_sample.F_MISS <= float(wildcards.st)][['SAMPLE']].to_csv(
            output.sample, index=False, header=False
        )


rule flt_VS:
    # filter variants and samples at the same time
    threads: 3
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{pfx}/vcfs/{reg}.vcf.gz.csi",
        variants = "{pfx}/variants_{vt}.txt",
        samples = "{pfx}/samples_{st}.txt",
    output:
        vcf = "{pfx}-V{vt}_S{st}/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -R {input.variants} {input.vcf} '
        '| bcftools view -S {input.samples} -o {output.vcf}'


rule prep_FWS:
    localrule: True
    input:
        fws = f'{{pfx}}/concat/fws-{complete_region}/fws.tsv'
    output:
        fws = '{pfx}/fws_{fws}.txt'
    run:
        import pandas as pd
        df_fws = pd.read_table(input.fws)
        df_fws[df_fws.FWS >= float(wildcards.fws)].SAMPLE.to_csv(
            output.fws, index=False, header=False,
        )


rule flt_FWS:
    # filter for Fws
    threads:2
    input:
        vcf = '{pfx}/vcfs/{reg}.vcf.gz',
        vcf_idx = '{pfx}/vcfs/{reg}.vcf.gz.csi',
        samples = '{pfx}/fws_{fws}.txt'
    output:
        vcf = '{pfx}-FWS{fws}/vcfs/{reg}.vcf.gz'
    shell:
        'bcftools view -S {input.samples} -o {output} {input.vcf}'


# -- analysis and calculation rules --


rule check_QC:
    threads: 1
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
        #"{pfx}/vcfs/{reg}.zarr"
    output:
        sample = "{pfx}/QC/{reg}-samples.tsv",
        variant = "{pfx}/QC/{reg}-variants.tsv"
    shell:
        'spcli vcf2qc --outsample {output.sample} --outvariant {output.variant} {input}'


rule check_single_QC:
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
        #"{pfx}/{reg}.zarr"
    output:
        sample = "{pfx}/{reg}.sample-qc.tsv",
        variant = "{pfx}/{reg}.variant-qc.tsv"
    log:
        "{pfx}/logs/check_single_QC-{reg}.log"
    benchmark:
        "{pfx}/logs/benchmark-check_single_QC-{reg}.txt"
    shell:
        'spcli vcf2qc --fraction --outsample {output.sample} --outvariant {output.variant} {input} 2> {log}'


rule plot_single_QC_png:
    threads: 1
    input:
        "{pfx}/{fname}.{var}-qc.tsv"
    output:
        "{pfx}/{fname}.{var}-qc.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        df = pd.read_table(input[0], sep='\t')
        # plot
        ax = sns.scatterplot(x=range(len(df)),
                             y=df.F_MISS.sort_values(),
                             s=3,
                             linewidth=0,
                             alpha=1.0)
        ax.set_ylabel(f'F_MISS')
        ax.set_xlabel(f'{wildcards.var} index (n={len(df)})')
        plt.savefig(output[0])
        plt.close()


use rule plot_single_QC_png as plot_single_QC_pdf with:
    output:
        "{pfx}/{reg}.{var}-qc.pdf"


rule plot_single_QC_meta_png:
    localrule: True
    input:
        tsv = "{pfx}/{reg}.{var}-qc.tsv",
    output:
        "{pfx}/{reg}.{var}-qc-{grp}.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt
        # from seqpy.core.bioio import tabutils

        df = pd.read_table(input[0], sep='\t')
        # df_meta = tabutils.join_metafile(df.SAMPLE, meta_file, data=df)
        df_meta = df.merge(pd.read_table(meta_file), how='left')
        df_meta.sort_values(by=[wildcards.grp, 'F_MISS'], inplace=True)
        # plot
        ax = sns.scatterplot(x=range(len(df_meta)),
                             y=df_meta.F_MISS,
                             hue=df_meta[wildcards.grp],
                             palette=sns.color_palette(color_palettes["xgfs_normal12"]),
                             s=4, linewidth=0, alpha=1.0)
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set_ylabel(f'F_MISS')
        ax.set_xlabel(f'{wildcards.var} index (n={len(df)})')
        plt.tight_layout()
        plt.savefig(output[0])
        plt.close()


rule gather_sample_QC_png:
    localrule: True
    input:
        expand('{{pfx}}/QC/{reg}-samples.tsv', reg=REGIONS)
    output:
        outfile = "{pfx}/QC/samples.tsv",
        outplot = "{pfx}/QC/samples.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        # gather
        dfs = [pd.read_table(infile, sep='\t') for infile in input]
        df = pd.concat(dfs).groupby('SAMPLE').agg('sum').reset_index()
        df['F_MISS'] = df.MISS / df.VARIANTS
        df.to_csv(output.outfile, sep='\t', index=False)

        # plot
        ax = sns.scatterplot(x=range(len(df)),
                             y=df.F_MISS.sort_values(),
                             s=3,
                             linewidth=0,
                             alpha=1.0)
        ax.set_ylabel('Variant F_MISS')
        ax.set_xlabel('Sample Index')
        plt.savefig(output.outplot)
        plt.close()


use rule gather_sample_QC_png as rule_gather_sample_QC_pdf with:
    output:
        outplot = "{pfx}/QC/samples.pdf"


rule gather_variant_QC:
    localrule: True
    input:
        expand('{{pfx}}/QC/{reg}-variants.tsv', reg=REGIONS)
    output:
        outfile = "{pfx}/QC/variants.tsv",
        outplot = "{pfx}/QC/variants.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        dfs = [pd.read_table(infile, sep='\t') for infile in input]
        df = pd.concat(dfs)
        df['F_MISS'] = df.MISS / df.SAMPLES
        df.to_csv(output.outfile, sep='\t', index=False)

        # plot
        ax = sns.scatterplot(x=range(len(df)),
                             y=df.F_MISS.sort_values(),
                             s=3,
                             linewidth=0,
                             alpha=1.0)
        ax.set_ylabel('Sample F_MISS')
        ax.set_xlabel('Variant Index')
        plt.savefig(output.outplot)
        plt.close()


use rule gather_variant_QC as gather_variant_QC_pdf with:
    output:
        outplot = "{pfx}/QC/variants.pdf"


rule calc_fws:
    # calculate FWS from GDS file using R moimix driven by gds2fws.py
    threads: 2
    input:
        "{pfx}/{complete_region}.gds"
    output:
        "{pfx}/fws-{complete_region}/fws.tsv"
    shell:
        "spcli $PYS/wgs/gds2fws.py -o {output} {input}"


rule plot_fws:
    threads: 1
    input:
        "{pfx}/fws-{complete_region}/fws.tsv"
    output:
        "{pfx}/fws-{complete_region}/fws-{grp}.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        fws_df = pd.read_table(input[0])

        ax = sns.catplot(
            data=fws_meta_df,
            x=wildcards.grp,
            order=sorted(fws_meta_df[wildcards.grp].unique()),
            y='FWS',
            palette=sns.color_palette(color_palettes["xgfs_normal12"]),
        )
        plt.xticks(rotation=60, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        plt.savefig(output[0])
        plt.close()


rule plot_fws_group:
    threads: 1
    input:
        "{pfx}/fws-{complete_region}/fws.tsv"
    output:
        "{pfx}/fws-{complete_region}/fws-{grp}.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        fws_df = pd.read_table(input[0])
        meta_df = pd.read_table(meta_file)
        fws_meta_df = fws_df.merge(meta_df)

        ax = sns.catplot(
            data=fws_meta_df,
            x=wildcards.grp,
            order=sorted(fws_meta_df[wildcards.grp].unique()),
            y='FWS',
            palette=sns.color_palette(color_palettes["xgfs_normal12"]),
        )
        plt.xticks(rotation=60, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        plt.savefig(output[0])
        plt.close()


ruleorder: plot_fws_group > plot_fws

rule chk_monoclonal_samples:
    localrule: True
    threads: 1
    input:
        "{pfx}/fws-{complete_region}/fws.tsv"
    output:
        "{pfx}/fws-{complete_region}/monoclonal-samples.txt"
    run:
        import pandas as pd

        df = pd.read_table(input[0])
        df_monoclonal = df[df.FWS > 0.95]
        df_monoclonal.SAMPLE.to_csv(output[0], header=False, index=False)

# -- NJ and genetic-distance stuff --

rule calc_gendist:
    # this rule generate proportional genetic distance matrix from zarr file
    threads: 4
    input:
        "{pfx}/{complete_region}.zarr"
    output:
        "{pfx}/gendist-{complete_region}/distance.tsv"
    shell:
        "spcli $PYS/wgs/zarr2dist.py --mindepth {mindepth} --allelenumber 2 -o {output} {input}"


rule plot_njtree:
    # this rule generate NJ tree plot from genetic distance
    # requires:
    # - metafile
    threads: 1
    input:
        "{pfx}/gendist-{complete_region}/distance.tsv"
    output:
        "{pfx}/gendist-{complete_region}/njtree.pdf"
    shell:
        "Rscript $PYS/plt/dist2treeplot.R -o {output} {input}"


rule prep_dist_anno:
    threads: 1
    input:
        distmat = "{pfx}/gendist-{complete_region}/distance.tsv",
        meta_file = meta_file,
    output:
        indv = "{pfx}/gendist-{complete_region}/dist-{grp}.indv.tsv",
        group = "{pfx}/gendist-{complete_region}/dist-{grp}.group.tsv"
    shell:
        "spcli $PYS/wgs/tab2anno.py --metafile {input.meta_file}:SAMPLE,{wildcards.grp} "
        "--outsample {output.indv} --outgroup {output.group} {input.distmat}"


rule plot_njtree_meta:
    threads: 1
    input:
        dist = "{pfx}/gendist-{complete_region}/distance.tsv",
        indv = "{pfx}/gendist-{complete_region}/dist-{grp}.indv.tsv",
        group = "{pfx}/gendist-{complete_region}/dist-{grp}.group.tsv"
    output:
        "{pfx}/gendist-{complete_region}/njtree-G@{grp}.pdf"
    shell:
        "Rscript $PYS/plt/dist2treeplot.R -c {input.indv} -l {input.group} -o {output} {input.dist}"


rule plot_njtree_meta_type:
    threads: 1
    input:
        dist = "{pfx}/gendist-{complete_region}/distance.tsv",
        indv = "{pfx}/gendist-{complete_region}/dist-{grp}.indv.tsv",
        group = "{pfx}/gendist-{complete_region}/dist-{grp}.group.tsv"
    output:
        "{pfx}/gendist-{complete_region}/njtree-G@{grp}-T@{type}.pdf"
    shell:
        "Rscript $PYS/plt/dist2treeplot.R -c {input.indv} -l {input.group} -o {output} "
        "-t {wildcards.type} -L {input.dist}"

ruleorder: plot_njtree_meta_type > plot_njtree_meta


# -- hmmIBD stuff --

rule select_sample_by_population:
    input:
        #"{pfx}/fws-{complete_region}/monoclonal-samples.txt"
        sample = '{pfx}/{complete_region}.sample-qc.tsv'
    output:
        sample = "{pfx}/hmmIBD-{complete_region}/{population}.{grp}.sample.txt"
    run:
        import pandas as pd

        df = pd.read_table(input.sample, header=None)
        df.rename(columns={0: 'SAMPLE'}, inplace=True)
        meta_df = df.merge(pd.read_table(meta_file))
        pops = wildcards.population.split('+')
        selected_df = meta_df.loc[meta_df[wildcards.grp].isin(pops)]
        selected_df.SAMPLE.to_csv(output.sample, header=False, index=False)


rule zarr2hmmibd:
    threads: 1
    input:
        zarr = "{pfx}/{complete_region}.zarr",
        sample = "{pfx}/hmmIBD-{complete_region}/{population}.sample.txt"
    output:
        '{pfx}/hmmIBD-{complete_region}/{population}.hmmIBD-input.tsv'
    shell:
        "spcli zarr2hmmibd "
        "  -s {input.sample} "
        "  -o {output} "
        "  -d {mindepth} "
        "  -t {chromtranslation_file} "
        "  {input.zarr}"


rule hmmibd:
    threads: 1
    input:
        '{pfx}/{population}.hmmIBD-input.tsv'
    output:
        '{pfx}/{population}.hmm_fract.txt'
    log:
        '{pfx}/logs/hmmIBD-{population}.log'
    shell:
        "hmmIBD -i {input} -m 1000 -n 1000 -o {wildcards.pfx}/{wildcards.population} 2> {log}"


rule ibd_kde:
    threads: 1
    input:
        '{pfx}/{population}.{grp}.hmm_fract.txt'
    output:
        '{pfx}/{population}.{grp}.kde.png'
    shell:
        "spcli $PYS/plt/kdeplot.py -o {output} -t {wildcards.population} {input}:fract_sites_IBD"


rule plot_ibd:
    threads: 1
    input:
        '{pfx}/{population}.hmm_fract.txt'
    output:
        plot = '{pfx}/{population}.ibd_plot.pdf',
        cluster = '{pfx}/{population}.clusters.yaml'
    shell:
        "spcli $PYS/plt/hmmIBD2ig.py -o {output.plot} --outcluster {output.cluster} "
        "{input}:sample1,sample2,fract_sites_IBD"


rule plot_ibd_meta:
    threads: 1
    input:
        '{pfx}/{population}.hmm_fract.txt'
    output:
        plot = '{pfx}/{population}.ibd_plot.{grp}.pdf',
        cluster = '{pfx}/{population}.clusters.{grp}.yaml'
    shell:
        "spcli $PYS/plt/hmmIBD2ig.py -o {output.plot} --outcluster {output.cluster} "
        "--metafile {meta_file}:{wildcards.grp} {input}:sample1,sample2,fract_sites_IBD"    


rule independent_samples:
    localrule: True
    input:
        cluster = '{pfx}/hmmIBD-{complete_region}/{population}.clusters.yaml',
        sample_qc = '{pfx}/{complete_region}.sample-qc.tsv'
    output:
        '{pfx}/hmmIBD-{complete_region}/{population}.independent_samples.txt'
    shell:
        "spcli $PYS/gen/get-independent-samples.py -o {output} --sample-qc {input.sample_qc} {input.cluster}"


# -- miscelaneous utilities --

rule concat_vcfs:
    # concatenate VCF files into a single VCF file
    # requires:
    #    REGIONS
    #    complete_region [the name to assign for the new VCF file]
    threads: 2
    input:
        expand("{{pfx}}/vcfs/{reg}.vcf.gz", reg=REGIONS)
    output:
        f"{{pfx}}/concat/{complete_region}.vcf.gz"
    shell:
        "bcftools concat -o {output} {input}"


rule vcf2zarr:
    # this rule convert a VCF file to a Zarr storage as Zarr storage is faster to read
    # requires:
    #   ploidy
    #   max_alt_alleles
    threads: 2
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        directory("{pfx}/{reg}.zarr")
    params:
        ploidy = config.get('ploidy', 2),
        max_alt_alleles = config.get('max_alt_alleles', 8)
    shell:
        'spcli vcf2zarr --ploidy {params.ploidy} --max_alt_alleles {params.max_alt_alleles} --outfile {output} {input}'


rule vcf2gds:
    # convert VCF to GDS using R SeqArray library driven by vcf2gds.py
    threads: 2
    input:
        f"{{pfx}}/{complete_region}.vcf.gz"
    output:
        f"{{pfx}}/{complete_region}.gds"
    shell:
        "spcli $PYS/wgs/vcf2gds.py -o {output} {input}"


rule zarr2tped:
    # convert Zarr to TPED
    threads: 2
    input:
        '{pfx}/{fname}.zarr'
    output:
        '{pfx}/{fname}.tped'
    shell:
        'spcli $PYS/wgs/zarr2tped.py -o {output} {input}'


rule index_csi:
    # indexing VCF with .csi
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.csi"
    shell:
        "bcftools index --csi {input}"


rule index_tbi:
    # indexing VCF with .tbi
    threads: 1
    input:
        "{pfx}/{reg}.vcf.gz"
    output:
        "{pfx}/{reg}.vcf.gz.tbi"
    shell:
        "bcftools index --tbi {input}"


#ruleorder: 


# EOF
