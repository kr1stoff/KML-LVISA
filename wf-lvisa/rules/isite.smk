rule map_umi_cover:
    input:
        umi_bam=rules.umi_filter.output,
        map_bam=rules.filter_bam.output.bam,
    output:
        umi_cov="is/{sample}.umi.coverage",
        all_cov="is/{sample}.all.coverage",
        umi_bed="is/{sample}.umi.bed",
    benchmark:
        ".log/is/{sample}.map_umi_cover.bm"
    log:
        ".log/is/{sample}.map_umi_cover.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        # umi bed 区域
        bedtools bamtobed -i {input.umi_bam} | bedtools merge > {output.umi_bed} 2> {log}
        # uniq umi bam 深度. coverage 解释: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html#default-behavior
        bedtools coverage -a {output.umi_bed} -b {input.umi_bam} > {output.umi_cov} 2>> {log}
        # 总 bam 深度
        bedtools coverage -a {output.umi_bed} -b {input.map_bam} > {output.all_cov} 2>> {log}
        """


rule isite_cover:
    input:
        umi_cov=rules.map_umi_cover.output.umi_cov,
        all_cov=rules.map_umi_cover.output.all_cov,
    output:
        "is/{sample}.is.coverage",
    benchmark:
        ".log/is/{sample}.isite_cover.bm"
    log:
        ".log/is/{sample}.isite_cover.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        # 过滤 1个loci < 2 UMI，support reads < 10
        python {config[my_scripts]}/filter_and_stat_is_by_coverage.py {input.umi_cov} {input.all_cov} {output} 2> {log}
        """
