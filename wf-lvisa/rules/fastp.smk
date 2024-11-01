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
        config["conda"]["basic"]
    params:
    # FIXME -l 最小长度要根据对应试剂盒 3LTR 长度来调整，目前是 82
        "-q 15 -u 40 -l 82 --adapter_sequence_r2 AAATGGTCTGAGGGATCTCTAGT --cut_right "
        "--cut_window_size 4 --cut_mean_quality 20 --correction --umi --umi_loc=read2 --umi_len=10 ",
    threads: config["threads"]["low"]
    shell:
        "fastp -w {threads} {params} -j {output.j} -h {output.h} -o {output.o} -O {output.O} -i {input[0]} -I {input[1]} &> {log}"


rule qc_stat:
    input:
        expand("qc/fastp/{sample}.json", sample=config["samples"]),
    output:
        "qc/fastp/fastp.stats.tsv",
    benchmark:
        ".log/fastp/qc_stat.bm"
    log:
        ".log/fastp/qc_stat.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/fastp_all_samples_qc.py {output} {input} &> {log}"
