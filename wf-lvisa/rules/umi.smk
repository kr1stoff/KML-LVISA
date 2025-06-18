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
        # todo [241231] 1.每个UMI最小支持reads数目为3；2.分类UMI最多允许1个碱基差异
        # [250206] 艾吉泰康建议不要设置 UMI support reads 阈值, 默认 1
        "-d 1 -s 3",
    threads: config["threads"]["low"]
    conda:
        config["conda"]["basic2"]
    shell:
        """
        gencore {params} -i {input} -o {output.bam} -r {config[database][hg19]} -h {output.html} -j {output.json} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


rule umi_bam_to_bed:
    input:
        rules.umi_bam.output.bam,
    output:
        temp("umi/{sample}.raw.bed"),
    log:
        ".log/umi/{sample}.umi_bam_to_bed.log",
    benchmark:
        ".log/umi/{sample}.umi_bam_to_bed.bm"
    conda:
        config["conda"]["basic2"]
    shell:
        "bedtools bamtobed -i {input} > {output}"


rule filter_umi_raw_bed:
    input:
        rules.umi_bam_to_bed.output,
    output:
        fltr="umi/{sample}.filter.bed",
        waste="umi/{sample}.waste.bed",
    log:
        ".log/umi/{sample}.filter_umi_raw_bed.log",
    benchmark:
        ".log/umi/{sample}.filter_umi_raw_bed.bm"
    conda:
        config["conda"]["python"]
    params:
        read_threas=1,
    shell:
        # * [250206] 针对性能验证 ZQX4 单条非特异 read 噪音解决办法
        """
        python {config[my_scripts]}/filter_umi_raw_bed.py {input} {output.fltr} {output.waste} {params.read_threas} 2>> {log}
        """


rule umi_strand_qname_bed:
    input:
        rules.filter_umi_raw_bed.output.fltr,
    output:
        strand="umi/{sample}.strand.bed",
        qname="umi/{sample}.qname.bed",
    log:
        ".log/umi/{sample}.umi_strand_qname_bed.log",
    benchmark:
        ".log/umi/{sample}.umi_strand_qname_bed.bm"
    conda:
        config["conda"]["basic2"]
    params:
        # * [250110] bedtools merge 参数
        # -s 正反向链不合并
        # -c 4 保留 rname 提取 UMI, -o collapse 折叠所有 rname
        # -c 6 保留正反向链信息, -o distinct 去重正反向链
        strand="-s -c 6 -o distinct",
        qname="-s -c 4 -o collapse",
    shell:
        """
        bedtools merge {params.strand} -i {input} > {output.strand}
        bedtools merge {params.qname} -i {input} > {output.qname}
        """


rule merge_umi_bed:
    input:
        strand=rules.umi_strand_qname_bed.output.strand,
        qname=rules.umi_strand_qname_bed.output.qname,
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
