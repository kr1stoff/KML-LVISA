# 综合大报告的内容都在这个 smk 更新
rule effect_stat:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
    output:
        text="bigrprt/effect.txt",
        fig="bigrprt/effect.png",
    log:
        ".log/bigrprt/effect_stat.log",
    benchmark:
        ".log/bigrprt/effect_stat.bm",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_effect_stats.py {output.text} {output.fig} {input} &> {log}"


rule repeat_stat:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
    output:
        text="bigrprt/repeat.txt",
        fig="bigrprt/repeat.png",
    log:
        ".log/bigrprt/repeat_stat.log",
    benchmark:
        ".log/bigrprt/repeat_stat.bm",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_repeat.stats.py {output.text} {output.fig} {input} &> {log}"


rule oncogene_stat:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
    output:
        text="bigrprt/oncogene.tsv",
    log:
        ".log/bigrprt/oncogene_stat.log",
    benchmark:
        ".log/bigrprt/oncogene_stat.bm",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/total_onco_stats.py {output.text} {input} &> {log}"
