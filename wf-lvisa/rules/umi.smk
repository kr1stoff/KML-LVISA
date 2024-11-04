rule umi_bam:
    input:
        rules.filter_bam.output.bam,
    output:
        bam="umi/{sample}.bam",
        json="umi/{sample}.umi.json",
        html="umi/{sample}.umi.html",
        stat="umi/{sample}.umi.bam.stat",
    benchmark:
        ".log/umi/{sample}.umi_bam.bm"
    log:
        ".log/umi/{sample}.umi_bam.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        # umi gencore 貌似只能去重不能统计
        gencore -i {input} -o {output.bam} -r {config[database][hg19]} -h {output.html} -j {output.json} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


rule umi_filter:
    input:
        rules.umi_bam.output.bam,
    output:
        "umi/{sample}.umi.filter.bam",
    benchmark:
        ".log/umi/{sample}.umi_filter.bm"
    log:
        ".log/umi/{sample}.umi_filter.log",
    conda:
        config["conda"]["python"]
    shell:
        """
        # 过滤小于3条 reads 的 UMI
        python {config[my_scripts]}/filter_umi_bam.py {input} {output} 2> {log}
        """
