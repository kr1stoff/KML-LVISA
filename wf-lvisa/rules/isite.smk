rule map_all_cover:
    input:
        bed=rules.merge_umi_bed.output,
        bam=rules.filter_bam.output.bam,
    output:
        "isite/{sample}.all.coverage",
    benchmark:
        ".log/isite/{sample}.map_all_cover.bm"
    log:
        ".log/isite/{sample}.map_all_cover.log",
    conda:
        config["conda"]["basic2"]
    shell:
        """
        # 总 bam 深度
        bedtools coverage -a {input.bed} -b {input.bam} > {output} 2>> {log}
        """


rule isite_cover:
    input:
        bed=rules.merge_umi_bed.output,
        all_cov=rules.map_all_cover.output,
    output:
        "isite/{sample}.is.coverage",
    benchmark:
        ".log/isite/{sample}.isite_cover.bm"
    log:
        ".log/isite/{sample}.isite_cover.log",
    conda:
        config["conda"]["python"]
    shell:
        """
        # 过滤 support reads < 10
        python {config[my_scripts]}/filter_and_stat_is_by_coverage.py {input.bed} {input.all_cov} {output} 2> {log}
        """
