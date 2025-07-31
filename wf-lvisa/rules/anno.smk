rule anno_transcript_exon:
    input:
        rules.table_annovar.output,
    output:
        "anno/{sample}.is.tx_exon",
    params:
        config["database"]["refseq_longest_tx"],
    benchmark:
        ".log/anno/{sample}.transcript_exon.bm"
    log:
        ".log/anno/{sample}.transcript_exon.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/annovar_tx_exon.py"


rule anno_interpro_domain:
    input:
        rules.table_annovar.output,
    output:
        "anno/{sample}.is.interpro_domain",
    benchmark:
        ".log/anno/{sample}.transcript_exon.bm"
    log:
        ".log/anno/{sample}.transcript_exon.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/annovar_interpro_domain.py"


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
        inbed="anno/{sample}.input.bed",
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
        cut -f-3 {input} > {output.inbed}
        # CpG 要去重, CpG(10kb) 之间有 overlap. 后面都做去重
        bedtools intersect -wb -a {output.inbed} -b {config[database][cpg1kb]} | cut -f1,2,7 > {output.cpg1kb} 2> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][cpg2d5kb]} | cut -f1,2,7 > {output.cpg2d5kb} 2> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][cpg5kb]} | cut -f1,2,7 > {output.cpg5kb} 2> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][cpg10kb]} | cut -f1,2,7 > {output.cpg10kb} 2> {log}
        # TSS
        bedtools intersect -wb -a {output.inbed} -b {config[database][tss1kb]} | cut -f1,2,7 > {output.tss1kb} 2>> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][tss2d5kb]} | cut -f1,2,7 > {output.tss2d5kb} 2>> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][tss5kb]} | cut -f1,2,7 > {output.tss5kb} 2>> {log}
        bedtools intersect -wb -a {output.inbed} -b {config[database][tss10kb]} | cut -f1,2,7 > {output.tss10kb} 2>> {log}
        # Repeat
        # repName, repClass, repFamily
        bedtools intersect -wb -a {output.inbed} -b {config[database][rmsk]} | cut -f1,2,7 > {output.repeat} 2>> {log}
        # 每 1M Window GC 含量
        bedtools intersect -wb -a {output.inbed} -b {config[database][gc1mb]} | cut -f1,2,7 > {output.gc1mb} 2>> {log}
        """
