rule generate_anno_control:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
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
        anno=rules.comb_anno.output,
        control=rules.generate_anno_control.output,
    output:
        "anno-qc/{sample}.tsv",
    benchmark:
        ".log/anno/{sample}.anno_qc.bm"
    log:
        ".log/anno/{sample}.anno_qc.log",
    conda:
        config["conda"]["python"]
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/annotate_qc.py {input.control} {input.anno} {output} 2> {log}"


rule concat_all_anno:
    input:
        expand("anno-qc/{sample}.tsv", sample=config["samples"]),
    output:
        "anno-qc/all.anno.qc.xlsx",
    benchmark:
        ".log/anno/concat_all_anno.bm"
    log:
        ".log/anno/concat_all_anno.log",
    conda:
        config["conda"]["basic"]
    shell:
        "csvtk -t csv2xlsx {input} -o {output} 2> {log}"
