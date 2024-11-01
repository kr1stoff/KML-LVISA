rule fastqc:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        directory("qc/fastqc/{sample}"),
    benchmark:
        ".log/fastqc/{sample}.bm"
    log:
        ".log/fastqc/{sample}.log",
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    shell:
        "mkdir {output} && fastqc {input} -o {output} -t {threads} --extract &> {log}"
