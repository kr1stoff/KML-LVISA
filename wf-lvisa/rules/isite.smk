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
        "bedtools coverage -a {input.bed} -b {input.bam} > {output} 2>> {log}"


rule rmdup_cover:
    input:
        bed=rules.merge_umi_bed.output,
        bam=rules.rmdup_bam.output.view,
    output:
        "isite/{sample}.rmdup.coverage",
    benchmark:
        ".log/isite/{sample}.rmdup_cover.bm"
    log:
        ".log/isite/{sample}.rmdup_cover.log",
    conda:
        config["conda"]["basic2"]
    shell:
        "bedtools coverage -a {input.bed} -b {input.bam} > {output} 2> {log}"


rule isite_cover:
    input:
        bed=rules.merge_umi_bed.output,
        all_cov=rules.map_all_cover.output,
        rmdup_cov=rules.rmdup_cover.output,
    output:
        "isite/{sample}.is.coverage",
    benchmark:
        ".log/isite/{sample}.isite_cover.bm"
    log:
        ".log/isite/{sample}.isite_cover.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/filter_and_stat_is_by_coverage.py"
