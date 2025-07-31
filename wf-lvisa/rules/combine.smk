# 合并表格
rule combine:
    input:
        rules.isite_cover.output,
        rules.anno_effect.output,
        rules.anno_oncokb.output,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg1kb,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg2d5kb,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg5kb,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg10kb,
        rules.bedtools_anno_cpg_tss_repeat.output.tss1kb,
        rules.bedtools_anno_cpg_tss_repeat.output.tss2d5kb,
        rules.bedtools_anno_cpg_tss_repeat.output.tss5kb,
        rules.bedtools_anno_cpg_tss_repeat.output.tss10kb,
        rules.bedtools_anno_cpg_tss_repeat.output.repeat,
        rules.bedtools_anno_cpg_tss_repeat.output.gc1mb,
        rules.anno_full_name.output,
        rules.anno_transcript_exon.output,
        rules.anno_interpro_domain.output,
    output:
        "combine/{sample}.tsv",
    benchmark:
        ".log/combine/{sample}.comb_anno.bm"
    log:
        ".log/combine/{sample}.comb_anno.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/combine_annotation.py {input} {output} 2> {log}"
