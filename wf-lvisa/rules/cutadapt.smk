# fastp 切不干净 umi_adapter, 加一步 cutadapt;
# UMI 前的一段接头 R1 用正向 TGGATAAAGTCGGA, R2 用反向互补 TCCGACTTTATCCA
rule cutadapt:
    input:
        rules.fastp_pe.output.o,
        rules.fastp_pe.output.O,
    output:
        "qc/cutadapt/{sample}.1.fastq.gz",
        "qc/cutadapt/{sample}.2.fastq.gz",
    benchmark:
        ".log/cutadapt/{sample}.cutadapt.bm"
    log:
        ".log/cutadapt/{sample}.cutadapt.log",
    conda:
        config["conda"]["basic"]
    params:
        # [20260204] cutadapt 参数
        # -G: R2 移除 5' 端; -a: R1 移除 3' 端; ^ 表示在 read 最头部; $ 表示在 read 最尾部;
        "-G ^TCCGACTTTATCCA -a TGGATAAAGTCGGA$",
    shell:
        """
        cutadapt {params} -o {output[0]} -p {output[1]} {input[0]} {input[1]} &> {log}
        """
