

# get params from config

known_variant_dir = config.get('known_variant_dir')
refseq = config.get('refseq', "")

REGIONS = config.get('regions', [])
complete_region = config.get('complete_region', 'complete')
# use empty list to indicate unused variables
variant_file = config.get('variant_file', [])
sample_file = config.get('sample_file', [])
meta_file = config.get('metadata_file', [])

#MAC = config.get('MAC', 3)
#source_dir = config.get('source_dir', None)
#sample_threshold = float(config.get('sample_threshold', 0.50))
#variant_threshold = float(config.get('variant_threshold', 0.50))
mindepth = int(config.get('mindepth', 5))


wildcard_constraints:
    reg = r'[\w-]+',
    MAC = r'\d+',
    mindepth = r'\d+',
    st = r'[\.\d]+',
    vt = r'[\.\d]+',
    var = r'[\w]+',


# -- sample and variant handling --


rule flt_pass:
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        '{pfx}-PASS-d{mindepth}/vcfs/{reg}.vcf.gz'
    shell:
        'bcftools view -R {known_variant_dir}/{wildcards.reg}.bed.gz {input} '
        '| bcftools +setGT -- -t q -n . -i "FORMAT/DP<{wildcards.mindepth}" '
        '| bcftools view -o {output}'


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
        #'bcftools +fill-tags {input} -- -t AF '
        #'| bcftools view -e "MAF<{wildcards.MAF}" -o {output} {input}'
        'bcftools +fill-tags {input} -- -t AF '
        '| bcftools view -q {wildcards.MAF}:minor  -o {output} {input}'


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
        "gunzip -c {input} | spcli $PYS/wgs/vcf2dedup.py --keep-max-mac - 2> {log} | bgzip -c > {output}"


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
    shell:
        'bcftools view -S {input.sample_file} -o {output} {input.vcf}'


rule flt_variants:
    # filter VCFs by variants
    threads: 2
    input:
        vcf = "{src_dir}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{src_dir}/vcfs/{reg}.vcf.gz",
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
        df_variant[df_variant.F_MISS <= float(wildcards.vt)][['CHROM', 'POS']].to_csv(output.variant, sep='\t', index=False, header=False)


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
        df_sample[df_sample.F_MISS <= float(wildcards.st)][['SAMPLE']].to_csv(output.sample, index=False, header=False)


rule flt_VS:
    # filter variants and samples at the same time
    threads: 3
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx = "{pfx}/vcfs/{reg}.vcf.gz.csi",
        variant = "{pfx}/variants_{vt}.txt",
        samples = "{pfx}/samples_{st}.txt",
    output:
        vcf = "{pfx}-V{vt},S{st}/vcfs/{reg}.vcf.gz"
    shell:
        'bcftools view -R {input.variant} {input.vcf} '
        '| bcftools view -S {input.sample} -o {output.vcf}'


# -- analysis and calculation rules --


rule check_QC:
    threads: 1
    input:
        "{pfx}/vcfs/{reg}.zarr"
    output:
        sample = "{pfx}/QC/{reg}-samples.tsv",
        variant = "{pfx}/QC/{reg}-variants.tsv"
    shell:
        'spcli zarr2qc --mindepth {mindepth} --outsample {output.sample} --outvariant {output.variant} {input}'


rule check_single_QC:
    threads: 1
    input:
        "{pfx}/{reg}.zarr"
    output:
        sample = "{pfx}/{reg}.sample-qc.tsv",
        variant = "{pfx}/{reg}.variant-qc.tsv"
    log:
        "{pfx}/logs/check_single_QC-{reg}.log"
    benchmark:
        "{pfx}/logs/benchmark-check_single_QC-{reg}.txt"
    shell:
        'spcli zarr2qc --fraction --mindepth {mindepth} --outsample {output.sample} --outvariant {output.variant} {input} 2> {log}'


rule plot_single_QC_png:
    localrule: True
    input:
        "{pfx}/{reg}.{var}-qc.tsv"
    output:
        "{pfx}/{reg}.{var}-qc.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        df = pd.read_table(input[0], sep='\t')
        # plot
        ax = sns.scatterplot(x=range(len(df)), y=df.F_MISS.sort_values(), s=3, linewidth=0, alpha=1.0)
        ax.set_ylabel(f'F_MISS')
        ax.set_xlabel(f'{wildcards.var} index (n={len(df)})')
        plt.savefig(output[0])


use rule plot_single_QC_png as plot_single_QC_pdf with:
    output:
        "{pfx}/{reg}.{var}-qc.pdf"


rule plot_single_QC_meta_png:
    localrule: True
    input:
        tsv = "{pfx}/{reg}.{var}-qc.tsv",
        meta = meta_file.split(':')[0] if meta_file else []
    output:
        "{pfx}/{reg}.{var}-qc-{grp}.png"
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt
        from seqpy.core.bioio import tabutils

        df = pd.read_table(input[0], sep='\t')
        df_meta = tabutils.join_metafile(df.SAMPLE, meta_file, data=df)
        df_meta.sort_values(by=[wildcards.grp, 'F_MISS'], inplace=True)
        # plot
        ax = sns.scatterplot(x=range(len(df_meta)),
                             y=df_meta.F_MISS,
                             hue=df_meta[wildcard.grp],
                             palette='tab20',
                             s=4, linewidth=0, alpha=1.0)
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set_ylabel(f'F_MISS')
        ax.set_xlabel(f'{wildcards.var} index (n={len(df)})')
        plt.tight_layout()
        plt.savefig(output[0])


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
        df = pd.concat(dfs).groupby('SAMPLE').agg(sum).reset_index()
        df['F_MISS'] = df.MISS / df.VARIANTS
        df.to_csv(output.outfile, sep='\t', index=False)

        # plot
        ax = sns.scatterplot(x=range(len(df)), y=df.F_MISS.sort_values(), s=3, linewidth=0, alpha=1.0)
        ax.set_ylabel('Variant F_MISS')
        ax.set_xlabel('Sample Index')
        plt.savefig(output.outplot)


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
        ax = sns.scatterplot(x=range(len(df)), y=df.F_MISS.sort_values(), s=3, linewidth=0, alpha=1.0)
        ax.set_ylabel('Sample F_MISS')
        ax.set_xlabel('Variant Index')
        plt.savefig(output.outplot)


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


# EOF
