rule diversity:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.is.diversity.tsv",
    benchmark:
        ".log/stats/{sample}.diversity.bm"
    log:
        ".log/stats/{sample}.diversity.log",
    conda:
        config["conda"]["R"]
    shell:
        "Rscript {config[my_scripts]}/is_diversity.R {input} {output} 2> {log}"


rule density_plot:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.chromosome_density.png",
    benchmark:
        ".log/stats/{sample}.density.bm"
    log:
        ".log/stats/{sample}.density.log",
    conda:
        config["conda"]["R"]
    shell:
        "Rscript {config[my_scripts]}/CMplot.R {input} {output} 2> {log}"


rule top10_is:
    input:
        rules.comb_anno.output,
    output:
        top10tsv="stats/{sample}.top10_is.tsv",
        top10png="stats/{sample}.top10_is.png",
    benchmark:
        ".log/stats/{sample}.top10_is.bm"
    log:
        ".log/stats/{sample}.top10_is.log",
    shell:
        "python {config[my_scripts]}/top10_is.py {input} {output.top10tsv} {output.top10png} 2> {log}"


rule chrom_dist:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.chrom_dist.png",
    benchmark:
        ".log/stats/{sample}.chrom_dist.bm"
    log:
        ".log/stats/{sample}.chrom_dist.log",
    shell:
        "python {config[my_scripts]}/chromosome_distribute_barplot.py {config[database][hg19comp]} {input} {output} 2> {log}"


rule effect_plot:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.effect_pie.png",
    benchmark:
        ".log/stats/{sample}.effect_plot.bm"
    log:
        ".log/stats/{sample}.effect_plot.log",
    shell:
        "python {config[my_scripts]}/effect_pie.py {input} {output} 2> {log}"


rule repclass_plot:
    input:
        rules.comb_anno.output,
    output:
        "stats/{sample}.repclass_pie.png",
    benchmark:
        ".log/stats/{sample}.repclass_plot.bm"
    log:
        ".log/stats/{sample}.repclass_plot.log",
    shell:
        "python {config[my_scripts]}/reclass_pie.py {input} {output} 2> {log}"
