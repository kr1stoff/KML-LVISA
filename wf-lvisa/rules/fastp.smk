rule fastp_pe:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        j="qc/fastp/{sample}.json",
        h="qc/fastp/{sample}.html",
        o="qc/fastp/{sample}.1.fastq.gz",
        O="qc/fastp/{sample}.2.fastq.gz",
    benchmark:
        ".log/fastp/{sample}.fastp_pe.bm"
    log:
        ".log/fastp/{sample}.fastp_pe.log",
    conda:
        config["conda"]["basic2"]
    params:
        # 艾吉泰康
        # 1. UMI 组合 384 个, 长度 8 bp
        # 2. 3LTR 长度 60 bp + 最小比对长度 20 bp
        # ! 3. [20260203] 切掉R2 UMI 拆分后的3'端 14 bp umi_dapter 序列
        # 如果测穿 R1 3'端是 adapter+UMI, R2 3'端是 LTR 序列, 在 cutadapt 处理
        "-q 15 -u 40 -l 80 --trim_front2 14 --umi --umi_loc=read2 --umi_len=8",
    threads: config["threads"]["low"]
    shell:
        """
        # UMI, R2测穿反向3LTR
        fastp -w {threads} {params} \
            -j {output.j} -h {output.h} \
            -o {output.o} -O {output.O} \
            -i {input[0]} -I {input[1]} &> {log}
        """


rule qc_stat:
    input:
        expand("qc/fastp/{sample}.json", sample=config["samples"]),
    output:
        "qc/fastp/fastp.stats.tsv",
        "qc/fastp/fastp.stats.xlsx",
    benchmark:
        ".log/fastp/qc_stat.bm"
    log:
        ".log/fastp/qc_stat.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/fastp_all_samples_qc.py {output} {input} &> {log}"
