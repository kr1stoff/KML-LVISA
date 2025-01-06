rule fq2fa:
    input:
        "qc/fastp/{sample}.1.fastq.gz",
    output:
        # temp() 临时文件, 会在流程结束后删除
        temp("3ltr/{sample}.fa"),
    benchmark:
        ".log/3ltr/{sample}.fq2fa.bm"
    log:
        ".log/3ltr/{sample}.fq2fa.log",
    conda:
        config["conda"]["basic"]
    shell:
        "seqtk seq -A {input} > {output} 2> {log}"


rule blastn_ltr:
    input:
        rules.fq2fa.output,
    output:
        temp("3ltr/{sample}.blast.out"),
    benchmark:
        ".log/3ltr/{sample}.blastn_ltr.bm"
    log:
        ".log/3ltr/{sample}.blastn_ltr.log",
    conda:
        config["conda"]["basic"]
    params:
        "-task blastn-short -word_size 7 -max_hsps 100 -max_target_seqs 1000 "
        "-outfmt '6 length nident pident sstrand sseqid sstart send slen qseqid qstart qend qlen qcovs evalue bitscore'",
    threads: config["threads"]["high"]
    shell:
        "blastn {params} -num_threads {threads} -query {input} -db {config[database][3ltr]} -out {output}  &> {log}"


rule get_header_from_blast:
    input:
        rules.blastn_ltr.output,
    output:
        temp("3ltr/{sample}.qseqid_read1.txt"),
    benchmark:
        ".log/3ltr/{sample}.get_header_from_blast.bm"
    log:
        ".log/3ltr/{sample}.get_header_from_blast.log",
    shell:
        "python {config[my_scripts]}/get_fastq_header_from_blast.py {input} {output} &> {log}"


rule grep_fq_by_header:
    input:
        rules.get_header_from_blast.output,
        "qc/fastp/{sample}.1.fastq.gz",
        "qc/fastp/{sample}.2.fastq.gz",
    output:
        temp("3ltr/{sample}.1.3ltr.fq"),
        temp("3ltr/{sample}.2.3ltr.fq"),
    benchmark:
        ".log/3ltr/{sample}.grep_fq_by_header.bm"
    log:
        ".log/3ltr/{sample}.grep_fq_by_header.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        seqkit grep -f {input[0]} {input[1]} -o {output[0]} &> {log}
        seqkit grep -f {input[0]} {input[1]} -o {output[1]} &>> {log}
        """
