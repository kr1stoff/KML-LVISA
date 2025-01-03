rule snpeff:
    input:
        rules.isite_cover.output,
    output:
        "anno/{sample}.is.bed.snpEff",
    benchmark:
        ".log/anno/{sample}.snpeff.bm"
    log:
        ".log/anno/{sample}.snpeff.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        snpEff -dataDir {config[database][snpeff]} -i bed -chr chr -geneId -canon -noStats hg19 {input} > {output} 2>> {log}
        """


rule anno_effect_oncokb:
    input:
        rules.snpeff.output,
    output:
        effect="anno/{sample}.is.effect",
        oncokb="anno/{sample}.is.oncokb",
    benchmark:
        ".log/anno/{sample}.anno_effect_oncokb.bm"
    log:
        ".log/anno/{sample}.anno_effect_oncokb.log",
    shell:
        """
        # 注释 position, gene, effect
        python {config[my_scripts]}/annotate_effect.py {input} {output.effect} 2> {log}
        # 注释 oncokb
        python {config[my_scripts]}/annotate_oncokb.py {output.effect} {output.oncokb} {config[database][oncokb]} 2>> {log}
        """


rule bedtools_anno_cpg_tss_repeat:
    input:
        rules.isite_cover.output,
    output:
        cpg="anno/{sample}.is.cpg",
        tss="anno/{sample}.is.tss",
        repeat="anno/{sample}.is.repeat",
    benchmark:
        ".log/anno/{sample}.bedtools_anno_cpg_tss_repeat.bm"
    log:
        ".log/anno/{sample}.bedtools_anno_cpg_tss_repeat.log",
    conda:
        config["conda"]["basic"]
    shell:
        """
        # CpG 要去重, CpG(10kb) 之间有 overlap. 后面都做去重
        bedtools intersect -wb -a {input} -b {config[database][cpg10kb]} | cut -f1,2,9 > {output.cpg} 2> {log}
        # TSS
        bedtools intersect -wb -a {input} -b {config[database][switchDbTss10kb]} | cut -f1,2,9 > {output.tss} 2>> {log}
        # Repeat
        # repName, repClass, repFamily
        bedtools intersect -wb -a {input} -b {config[database][rmsk]} | cut -f1,2,9 > {output.repeat} 2>> {log}
        """


rule comb_anno:
    input:
        rules.isite_cover.output,
        rules.anno_effect_oncokb.output.effect,
        rules.anno_effect_oncokb.output.oncokb,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg,
        rules.bedtools_anno_cpg_tss_repeat.output.tss,
        rules.bedtools_anno_cpg_tss_repeat.output.repeat,
    output:
        "anno/{sample}.is.combine.tsv",
    benchmark:
        ".log/anno/{sample}.comb_anno.bm"
    log:
        ".log/anno/{sample}.comb_anno.log",
    shell:
        """
        # 合并注释结果
        python {config[my_scripts]}/combine_annotation.py {input} {output} 2> {log}
        """
