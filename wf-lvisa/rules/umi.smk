rule umi_bam:
    input:
        rules.filter_bam.output.bam,
    output:
        bam="umi/{sample}.bam",
        sam="umi/{sample}.sam",
        json="umi/{sample}.umi.json",
        html="umi/{sample}.umi.html",
        stat="umi/{sample}.umi.bam.stat",
    benchmark:
        ".log/umi/{sample}.umi_bam.bm"
    log:
        ".log/umi/{sample}.umi_bam.log",
    params:
        # * [241231] 1.每个UMI最小支持reads数目为3；2.分类UMI最多允许1个碱基差异
        "-s 3 -d 1",
    threads: config["threads"]["low"]
    conda:
        config["conda"]["basic"]
    shell:
        """
        gencore {params} -i {input} -o {output.bam} -r {config[database][hg19]} -h {output.html} -j {output.json} 2> {log}
        samtools view -@ {threads} {output.bam} > {output.sam}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


rule umi_filter:
    input:
        sam=rules.umi_bam.output.sam,
        umi=config["database"]["umi"],
    output:
        "umi/{sample}.umi.bed",
    benchmark:
        ".log/umi/{sample}.umi_filter.bm"
    log:
        ".log/umi/{sample}.umi_filter.log",
    conda:
        config["conda"]["python"]
    shell:
        """
        python {config[my_scripts]}/filter_umi_bam.py {input.sam} {input.umi} {output} 2>> {log}
        """
