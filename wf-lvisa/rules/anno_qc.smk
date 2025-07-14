rule generate_anno_control:
    input:
        expand("anno/{sample}.is.combine.tsv", sample=config["samples"]),
    output:
        "anno-qc/controls.tsv",
    params:
        pc_site_file=config["database"]["pc"],
    benchmark:
        ".log/anno-qc/generate_anno_control.bm"
    log:
        ".log/anno-qc/generate_anno_control.log",
    conda:
        config["conda"]["python"]
    shell:
        "python {config[my_scripts]}/generate_anno_control.py {output} {params.pc_site_file} {input} 2> {log}"


rule anno_qc:
    input:
        anno=rules.comb_anno.output,
        control=rules.generate_anno_control.output,
    output:
        "anno-qc/{sample}.anno.qc.tsv",
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
