

# -- signal of selection scan soss_ --

rule soss_prep:
    localrule: True
    # prepare the file so the steps can be run locally within selscan-{complete_region} directory
    input:
        '{pfx}/hmmIBD-{complete_region}/{population}.independent_samples.txt'
    output:
        '{pfx}/selscan-{complete_region}/{population}.sample-list.txt'
    shell:
        "ln -sr {input} {output}"

rule soss_prep_vcf:
    threads: 2
    input:
        vcf = '{pfx}/{complete_region}.vcf.gz',
        idx = '{pfx}/{complete_region}.vcf.gz.tbi',
        sample = '{pfx}/selscan-{complete_region}/{population}.sample-list.txt'
    output:
        vcf = '{pfx}/selscan-{complete_region}/{population}.vcf.gz'
    shell:
        'bcftools view -S {input.sample} {input.vcf} | '
        'spcli $PYS/wgs/vcf2sethet.py -o {output.vcf} '
        '  --min-alt-count 2 --min-alt-ratio 0.45 --set-het-to-ref --set-missing-to-ref'

rule index_vcf_region:
    threads: 1
    input:
        vcf = '{pfx}/selscan-{complete_region}/{population}.vcf.gz',
    output:
        idx = '{pfx}/selscan-{complete_region}/{population}.vcf.gz.tbi'
    shell:
        'bcftools index --tbi {input.vcf}'

rule soss_prep_vcf_region:
    localrule: True
    #threads: 1
    input:
        vcf = '{pfx}/selscan-{complete_region}/{population}.vcf.gz',
        idx = '{pfx}/selscan-{complete_region}/{population}.vcf.gz.tbi',
    output:
        vcf = '{pfx}/selscan-{complete_region}/{pop1}-+-{pop2}-xpnsl/{population}-{reg}.vcf.gz'
    shell:
        'bcftools view -o {output.vcf} -r {wildcards.reg} {input.vcf}'


rule soss_vcf2tped:
    # convert vcf to tpepd based on chromosomes
    threads: 1
    input:
        vcf = '{pfx}/{population}.vcf.gz'
    output:
        tped = '{pfx}/{pop1}-+-{pop2}-xpehh/{population}-{reg}.tped'
    shell:
        'spcli $PYS/wgs/vcf2tped.py '
        '  --translation-file {chromtranslation_file} '
        '  --region {wildcards.reg} '
        '  -o {output.tped} '
        '  {input.vcf}'


rule soss_selscan_xpehh:
    # options to use --wagh
    threads: 4
    input:
        pop1 = '{pfx}/{pop1}-+-{pop2}-xpehh/{pop1}-{reg}.tped',
        pop2 = '{pfx}/{pop1}-+-{pop2}-xpehh/{pop2}-{reg}.tped',
    output:
        outfile = '{pfx}/{pop1}-+-{pop2}-xpehh/{reg}.xpehh.out',
    shell:
        'selscan --xpehh  --threads {threads} '
        '  --out {wildcards.pfx}/{wildcards.pop1}-+-{wildcards.pop2}-xpehh/{wildcards.reg} '
        '  --tped {input.pop1} --tped-ref {input.pop2}'


rule soss_concat_normalize_xpehh:
    localrule: True
    input:
        expand('{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{reg}.xpehh.out', reg=REGIONS)
    output:
        outfile = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{complete_region}.xpehh.tsv'
    params:
        output = lambda w, input: ' '.join(f'{infile}.norm' for infile in input)
    shell:
        'norm --xpehh --files {input} && '
        'cat {params.output} > {output}'


rule soss_selscan_xpnsl:
    threads: 4
    input:
        pop1 = '{pfx}/{pop1}-+-{pop2}-xpnsl/{pop1}-{reg}.vcf.gz',
        pop2 = '{pfx}/{pop1}-+-{pop2}-xpnsl/{pop2}-{reg}.vcf.gz',
    output:
        outfile = '{pfx}/{pop1}-+-{pop2}-xpnsl/{reg}.xpnsl.out',
    shell:
        'selscan --xpnsl  --threads {threads} --unphased'
        '  --out {wildcards.pfx}/{wildcards.pop1}-+-{wildcards.pop2}-xpnsl/{wildcards.reg} '
        '  --vcf {input.pop1} --vcf-ref {input.pop2}'


rule soss_concat_normalize_xpnsl:
    localrule: True
    input:
        expand('{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{reg}.xpnsl.out', reg=REGIONS)
    output:
        outfile = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{complete_region}.xpnsl.tsv'
    params:
        output = lambda w, input: ' '.join(f'{infile}.norm' for infile in input)
    shell:
        'norm --xpnsl --files {input} && '
        'cat {params.output} > {output}'


rule soss_plot_xpnsl:
    localrule: True
    input:
        xpnsl =  f'{{pfx}}/{complete_region}.xpnsl.tsv'
    output:
        plot = f'{{pfx}}/{complete_region}.xpnsl.png'
    shell:
        'spcli $PYS/plt/mhtplot2.py -o {output.plot} --column normxpnsl --use-id id'

# EOF
