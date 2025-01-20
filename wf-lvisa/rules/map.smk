rule map2hg19:
    input:
        rules.grep_fq_by_header.output,
    output:
        bam="map/{sample}.bam",
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
        bwa mem -t {threads} {params.bwa} {config[database][hg19]} {input} | \
            samtools view -@ {threads} {params.view} - | \
            samtools sort -@ {threads} -o {output.bam} - 2> {log}
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
    params:
        # ! 过滤
        # 1.mapq < 30
        # 2.map length < 55
        # 3. flag 4 8 256 512 2048
        # 4. TLEN > 1000
        # 5. read_length (CIGAR M or S) > 50
        "-hbS -q 30 -m 55 -F 2828 -e 'tlen < 1000 && qlen-sclen > 50'",
    threads: config["threads"]["high"]
    shell:
        """
        samtools view -@ {threads} {params} {input} -o {output.bam} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """
