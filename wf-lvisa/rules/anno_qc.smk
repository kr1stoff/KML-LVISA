# 新增阳控/阴控/批次质控
rule generate_anno_control:
    input:
        expand("combine/{sample}.tsv", sample=config["samples"]),
    output:
        "anno-qc/controls.tsv",
    benchmark:
        ".log/anno-qc/generate_anno_control.bm"
    log:
        ".log/anno-qc/generate_anno_control.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/generate_anno_control.py {output} {input} 2> {log}"


rule anno_qc:
    input:
        anno=rules.combine.output,
        control=rules.generate_anno_control.output,
    output:
        "anno-qc/{sample}.qc.tsv",
    benchmark:
        ".log/anno/{sample}.anno_qc.bm"
    log:
        ".log/anno/{sample}.anno_qc.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/annotate_qc.py {input.control} {input.anno} {output} 2> {log}"


rule anno_filter:
    input:
        rules.anno_qc.output,
    output:
        "anno-qc/{sample}.filter.tsv",
    benchmark:
        ".log/anno/{sample}.anno_filter.bm"
    log:
        ".log/anno/{sample}.anno_filter.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/anno_filter.py"


rule concat_all_anno:
    input:
        expand("anno-qc/{sample}.filter.tsv", sample=config["samples"]),
    output:
        "anno-qc/all.anno.xlsx",
    benchmark:
        ".log/anno/concat_all_anno.bm"
    log:
        ".log/anno/concat_all_anno.log",
    conda:
        config["conda"]["basic"]
    shell:
        "csvtk csv2xlsx -t -f {input} -o {output} 2> {log}"
