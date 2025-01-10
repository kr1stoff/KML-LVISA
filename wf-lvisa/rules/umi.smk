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
    params:
        # * [241231] 1.每个UMI最小支持reads数目为3；2.分类UMI最多允许1个碱基差异
        "-s 3 -d 1",
    threads: config["threads"]["low"]
    conda:
        config["conda"]["basic"]
    shell:
        """
        gencore {params} -i {input} -o {output.bam} -r {config[database][hg19]} -h {output.html} -j {output.json} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


rule umi_bam_to_bed:
    input:
        rules.umi_bam.output.bam,
    output:
        strand="umi/{sample}.strand.bed",
        qname="umi/{sample}.qname.bed",
    log:
        ".log/umi/{sample}.umi_bam_to_bed.log",
    benchmark:
        ".log/umi/{sample}.umi_bam_to_bed.bm"
    conda:
        config["conda"]["basic"]
    params:
        # * [250110] bedtools merge 参数
        # -s 正反向链不合并
        # -c 4 保留 rname 提取 UMI, -o collapse 折叠所有 rname
        # -c 6 保留正反向链信息, -o distinct 去重正反向链
        strand="-s -c 6 -o distinct",
        qname="-s -c 4 -o collapse",
    shell:
        """
        bedtools bamtobed -i {input} | bedtools merge {params.strand} > {output.strand}
        bedtools bamtobed -i {input} | bedtools merge {params.qname} > {output.qname}
        """


rule merge_umi_bed:
    input:
        strand=rules.umi_bam_to_bed.output.strand,
        qname=rules.umi_bam_to_bed.output.qname,
    output:
        "umi/{sample}.umi.bed",
    benchmark:
        ".log/umi/{sample}.merge_umi_bed.bm"
    log:
        ".log/umi/{sample}.merge_umi_bed.log",
    conda:
        config["conda"]["python"]
    params:
        umi=config["database"]["umi"],
    shell:
        """
        python {config[my_scripts]}/merge_umi_bed.py {input.strand} {input.qname} {params.umi} {output} 2>> {log}
        """


# rule umi_filter:
#     input:
#         sam=rules.umi_bam.output.sam,
#         umi=config["database"]["umi"],
#     output:
#         "umi/{sample}.umi.bed",
#     benchmark:
#         ".log/umi/{sample}.umi_filter.bm"
#     log:
#         ".log/umi/{sample}.umi_filter.log",
#     conda:
#         config["conda"]["python"]
#     shell:
#         """
#         python {config[my_scripts]}/filter_umi_bam.py {input.sam} {input.umi} {output} 2>> {log}
#         """
