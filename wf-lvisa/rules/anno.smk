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
        "snpEff -dataDir {config[database][snpeff]} -i bed -chr chr -geneId -canon -noStats hg19 {input} > {output} 2>> {log}"


rule anno_effect:
    input:
        rules.snpeff.output,
    output:
        "anno/{sample}.is.effect",
    benchmark:
        ".log/anno/{sample}.anno_effect.bm"
    log:
        ".log/anno/{sample}.anno_effect.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/annotate_effect.py {input} {output} 2> {log}"


rule anno_oncokb:
    input:
        rules.anno_effect.output,
    output:
        "anno/{sample}.is.oncokb",
    benchmark:
        ".log/anno/{sample}.anno_oncokb.bm"
    log:
        ".log/anno/{sample}.anno_oncokb.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/annotate_oncokb.py {input} {output} {config[database][oncokb]} 2> {log}"


rule anno_full_name:
    input:
        rules.anno_effect.output,
    output:
        "anno/{sample}.is.fullname",
    benchmark:
        ".log/anno/{sample}.anno_full_name.bm"
    log:
        ".log/anno/{sample}.anno_full_name.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/annotate_fullname.py {input} {output} {config[database][hgnc]} 2> {log}"


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
        config["conda"]["basic2"]
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
        rules.anno_effect.output,
        rules.anno_oncokb.output,
        rules.bedtools_anno_cpg_tss_repeat.output.cpg,
        rules.bedtools_anno_cpg_tss_repeat.output.tss,
        rules.bedtools_anno_cpg_tss_repeat.output.repeat,
        rules.anno_full_name.output,
    output:
        "anno/{sample}.is.combine.tsv",
    benchmark:
        ".log/anno/{sample}.comb_anno.bm"
    log:
        ".log/anno/{sample}.comb_anno.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/combine_annotation.py {input} {output} 2> {log}"
