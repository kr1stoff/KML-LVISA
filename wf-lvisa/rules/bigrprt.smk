# 综合大报告的内容都在这个 smk 更新
rule effect_stat:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
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
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
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
        "python {config[my_scripts]}/total_repeat.stats.py {output.text} {output.fig} {output.summary} {input} &> {log}"


rule oncogene_stat:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
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
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
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
