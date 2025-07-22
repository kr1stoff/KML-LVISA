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
        rules.table_annovar.output,
    output:
        "anno/{sample}.is.effect",
    benchmark:
        ".log/anno/{sample}.anno_effect.bm"
    log:
        ".log/anno/{sample}.anno_effect.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/annotate_effect.py"


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
        cpg1kb="anno/{sample}.is.cpg1kb",
        cpg2d5kb="anno/{sample}.is.cpg2d5kb",
        cpg5kb="anno/{sample}.is.cpg5kb",
        cpg10kb="anno/{sample}.is.cpg10kb",
        tss1kb="anno/{sample}.is.tss1kb",
        tss2d5kb="anno/{sample}.is.tss2d5kb",
        tss5kb="anno/{sample}.is.tss5kb",
        tss10kb="anno/{sample}.is.tss10kb",
        repeat="anno/{sample}.is.repeat",
        gc1mb="anno/{sample}.is.gc",
    benchmark:
        ".log/anno/{sample}.bedtools_anno_cpg_tss_repeat.bm"
    log:
        ".log/anno/{sample}.bedtools_anno_cpg_tss_repeat.log",
    conda:
        config["conda"]["basic2"]
    shell:
        """
        # CpG 要去重, CpG(10kb) 之间有 overlap. 后面都做去重
        bedtools intersect -wb -a {input} -b {config[database][cpg1kb]} | cut -f1,2,9 > {output.cpg1kb} 2> {log}
        bedtools intersect -wb -a {input} -b {config[database][cpg2d5kb]} | cut -f1,2,9 > {output.cpg2d5kb} 2> {log}
        bedtools intersect -wb -a {input} -b {config[database][cpg5kb]} | cut -f1,2,9 > {output.cpg5kb} 2> {log}
        bedtools intersect -wb -a {input} -b {config[database][cpg10kb]} | cut -f1,2,9 > {output.cpg10kb} 2> {log}
        # TSS
        bedtools intersect -wb -a {input} -b {config[database][tss1kb]} | cut -f1,2,9 > {output.tss1kb} 2>> {log}
        bedtools intersect -wb -a {input} -b {config[database][tss2d5kb]} | cut -f1,2,9 > {output.tss2d5kb} 2>> {log}
        bedtools intersect -wb -a {input} -b {config[database][tss5kb]} | cut -f1,2,9 > {output.tss5kb} 2>> {log}
        bedtools intersect -wb -a {input} -b {config[database][tss10kb]} | cut -f1,2,9 > {output.tss10kb} 2>> {log}
        # Repeat
        # repName, repClass, repFamily
        bedtools intersect -wb -a {input} -b {config[database][rmsk]} | cut -f1,2,9 > {output.repeat} 2>> {log}
        # 每 1M Window GC 含量
        bedtools intersect -wb -a {input} -b {config[database][gc1mb]} | cut -f1,2,9 > {output.gc1mb} 2>> {log}
        """


rule comb_anno:
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
