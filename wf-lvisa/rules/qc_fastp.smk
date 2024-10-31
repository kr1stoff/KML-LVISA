rule fastp_pe:
    input:
        '.rawdata/{sample}_1.fastq.gz',
        '.rawdata/{sample}_2.fastq.gz',
    output:
        j='qc/fastp/{sample}.json',
        h='qc/fastp/{sample}.html',
        o='qc/fastp/{sample}.1.fastq.gz',
        O='qc/fastp/{sample}.2.fastq.gz'
    benchmark:
        '.log/fastp/{sample}.fastp_pe.bm'
    log:
        '.log/fastp/{sample}.fastp_pe.log'
    conda:
        config['conda']['basic']
    params:
        '-q 15 -u 40 -t 0 -G -n 5 -l 15 -y'
    threads:
        config['threads']['low']
    shell:
        'fastp -w {threads} {params} -j {output.j} -h {output.h} -o {output.o} -O {output.O} -i {input[0]} -I {input[1]} &> {log}'


rule qc_stat:
    input:
        expand('qc/fastp/{sample}.json',sample=config["samples"])
    output:
        'qc/fastp/fastp.stats.xls'
    benchmark:
        '.log/fastp/qc_stat.bm'
    log:
        '.log/fastp/qc_stat.log'
    conda:
        config['conda']['python']
    script:
        '../scripts/fastp_all_samples_qc.py'
