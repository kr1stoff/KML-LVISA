rule diversity_index:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.is.diversity.tsv",
    benchmark:
        ".log/stats/{sample}.diversity_index.bm"
    log:
        ".log/stats/{sample}.diversity_index.log",
    conda:
        config["conda"]["R"]
    shell:
        """
        Rscript {config[my_scripts]}/is_diversity.R {input} {output} 2> {log}
        """


rule top_chrom_effect_repclass:
    input:
        rules.comb_anno.output,
    output:
        top10tsv="stats/{sample}.top10_is.tsv",
        top10png="stats/{sample}.top10_is.png",
        chrom_dist="stats/{sample}.chrom_dist.png",
        effect="stats/{sample}.effect_pie.png",
        repclass="stats/{sample}.repclass_pie.png",
    benchmark:
        ".log/stats/{sample}.top_chrom_effect_repclass.bm"
    log:
        ".log/stats/{sample}.top_chrom_effect_repclass.log",
    conda:
        config["conda"]["python"]
    shell:
        """
        # Top 10 IS 位点柱状图
        python {config[my_scripts]}/top10_is.py {input} {output.top10tsv} {output.top10png} 2> {log}
        # 染色体插入位点占比柱状图 (均一化化染色体长度)
        python {config[my_scripts]}/chromosome_distribute_barplot.py {config[database][hg19comp]} {input} {output.chrom_dist} 2>> {log}
        # snpEff effect 饼图
        python {config[my_scripts]}/effect_pie.py {input} {output.effect} 2>> {log}
        # Repeat repClass 饼图
        python {config[my_scripts]}/reclass_pie.py {input} {output.repclass} 2>> {log}
        """
