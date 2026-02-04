# fastp 切不干净, 加一步 cutadapt; UMI前的一段接头R1用正向 TGGATAAAGTCGGA, R2 用反向互补 TCCGACTTTATCCA
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
        "-G TCCGACTTTATCCA -g TGGATAAAGTCGGA",
    threads: config["threads"]["low"]
    shell:
        """
        cutadapt {params} -o {output[0]} -p {output[1]} {input[0]} {input[1]} &> {log}
        """
