# 综合大报告的内容都在这个 smk 更新
rule effect_stat:
    input:
        expand("combine/{sample}.tsv", sample=config["samples"]),
    output:
        text="bigrprt/effect.txt",
        fig="bigrprt/effect.png",
        summary="bigrprt/effect_summary.tsv",
    log:
        ".log/bigrprt/effect_stat.log",
    benchmark:
        ".log/bigrprt/effect_stat.bm"
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_effect_stats.py {output.text} {output.fig} {output.summary} {input} &> {log}"


rule repeat_stat:
    input:
        expand("combine/{sample}.tsv", sample=config["samples"]),
    output:
        text="bigrprt/repeat.txt",
        fig="bigrprt/repeat.png",
        summary="bigrprt/repeat_summary.tsv",
    log:
        ".log/bigrprt/repeat_stat.log",
    benchmark:
        ".log/bigrprt/repeat_stat.bm"
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_repeat_stats.py {output.text} {output.fig} {output.summary} {input} &> {log}"


# * oncogene 部分保存解读参数过滤后的结果, 其他都用过滤前的数据
rule oncogene_stat:
    input:
        expand("anno-qc/{sample}.filter.tsv", sample=config["samples"]),
    output:
        text="bigrprt/oncogene.tsv",
    log:
        ".log/bigrprt/oncogene_stat.log",
    benchmark:
        ".log/bigrprt/oncogene_stat.bm"
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_onco_stats.py {output.text} {input} &> {log}"


rule summary_cpg_tss:
    input:
        expand("combine/{sample}.tsv", sample=config["samples"]),
    output:
        cpg="bigrprt/cpg_summary.tsv",
        tss="bigrprt/tss_summary.tsv",
    log:
        ".log/bigrprt/summary_cpg_tss.log",
    benchmark:
        ".log/bigrprt/summary_cpg_tss.bm"
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/summary_cpg_tss.py {output.cpg} {output.tss} {input} &> {log}"


rule summary_chromosome:
    input:
        expand("combine/{sample}.tsv", sample=config["samples"]),
    output:
        "bigrprt/chromosome_summary.tsv",
    log:
        ".log/bigrprt/summary_chromosome.log",
    benchmark:
        ".log/bigrprt/summary_chromosome.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_chromosome.py"


rule summary_diversity:
    input:
        expand("stats/{sample}.is.diversity.tsv", sample=config["samples"]),
    output:
        "bigrprt/diversity_summary.tsv",
    log:
        ".log/bigrprt/summary_diversity.log",
    benchmark:
        ".log/bigrprt/summary_diversity.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_diversity.py"
