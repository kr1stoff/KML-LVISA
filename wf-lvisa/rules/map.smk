rule map2hg19:
    input:
        rules.grep_fq_by_header.output,
    output:
        bam=temp("map/{sample}.bam"),
        stat="map/{sample}.bam.stat",
    benchmark:
        ".log/map/{sample}.map2hg19.bm"
    log:
        ".log/map/{sample}.map2hg19.log",
    conda:
        config["conda"]["basic"]
    params:
        bwa="-M -Y -R '@RG\\tID:{sample}\\tSM:{sample}'",
        view="-hbS",
    threads: config["threads"]["high"]
    shell:
        """
        bwa mem -t {threads} {params.bwa} {config[database][hg19]} {input} 2> {log} | \
            samtools view -@ {threads} {params.view} - 2>> {log} | \
            samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


rule filter_bam:
    input:
        rules.map2hg19.output.bam,
    output:
        bam="map/{sample}.filter.bam",
        stat="map/{sample}.filter.bam.stat",
    benchmark:
        ".log/map/{sample}.filter_bam.bm"
    log:
        ".log/map/{sample}.filter_bam.log",
    conda:
        config["conda"]["basic"]
    # ! 过滤
    # 1.mapq < 20
    # 2.min length < 55
    # 3.flag 2828:
    #   4(read unmapped)
    #   8(mate unmapped)
    #   256(not primary alignment)
    #   512(read fails platform/vendor quality checks)
    #   2048(supplementary alignment).
    # 4.[250711] TLEN > 600 (原参数 1000)
    # 5.alignment length < 55, (query length - softclip length)
    params:
        "-hbS -q 20 -m 55 -F 2828 -e 'tlen < 600 && qlen-sclen > 55'",
    threads: config["threads"]["high"]
    shell:
        """
        samtools view -@ {threads} {params} {input} -o {output.bam} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


# [20250730 MXFA] 新增去重 reads 数统计. 仅统计去重 reads 使用, 不在 umi 流程中
# ! [20251223] 先过滤后去重, reads 不成对, 没有 ms(mate score) tag, samtools markdup 会报错. 调整为先去重后过滤(该步骤与仅影响去重reads计数)
rule rmdup_bam:
    input:
        rules.map2hg19.output.bam,
    output:
        sort_name=temp("map/{sample}.sorted_by_name.bam"),
        fixmate=temp("map/{sample}.fixmate.bam"),
        sort_coodi=temp("map/{sample}.sorted_by_coodinate.bam"),
        rmdup=temp("map/{sample}.rmdup.bam"),
        stats="map/{sample}.rmdup_stats.txt",
        view=temp("map/{sample}.rmdup.filter.bam"),
    benchmark:
        ".log/map/{sample}.rmdup_bam.bm"
    log:
        ".log/map/{sample}.rmdup_bam.log",
    conda:
        config["conda"]["basic"]
    params:
        sortn="-n",
        fixmate="-r -c -m",
        markdup="-r -s",
        view="-hbS -q 20 -m 55 -F 2828 -e 'tlen < 600 && qlen-sclen > 55'",
    threads: config["threads"]["low"]
    shell:
        """
        samtools sort -@ 32 {params.sortn} {input} -o {output.sort_name} 2> {log}
        samtools fixmate -@ 32 {params.fixmate} {output.sort_name} {output.fixmate} 2>> {log}
        samtools sort -@ 32 {output.fixmate} -o {output.sort_coodi} 2>> {log}
        samtools markdup -@ 32 {params.markdup} -f {output.stats} {output.sort_coodi} {output.rmdup} 2>> {log}
        samtools view -@ {threads} {params.view} {output.rmdup} -o {output.view} 2>> {log}
        """
