rule prepare_annovar_input1:
    input:
        rules.isite_cover.output,
    output:
        "annovar/{sample}.bed.fa",
    benchmark:
        ".log/annovar/{sample}.prepare_annovar_input1.bm"
    log:
        ".log/annovar/{sample}.prepare_annovar_input1.log",
    conda:
        config["conda"]["basic2"]
    shell:
        "bedtools getfasta -fi {config[database][hg19]} -bed {input} -fo {output} 2> {log}"


rule prepare_annovar_input2:
    input:
        rules.prepare_annovar_input1.output,
    output:
        "annovar/{sample}.vcf",
    benchmark:
        ".log/annovar/{sample}.prepare_annovar_input2.bm"
    log:
        ".log/annovar/{sample}.prepare_annovar_input2.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/prepare_annovar_fa2vcf.py"


rule convert2annovar:
    input:
        rules.prepare_annovar_input2.output,
    output:
        "annovar/{sample}.avinput",
    benchmark:
        ".log/annovar/{sample}.convert2annovar.bm"
    log:
        ".log/annovar/{sample}.convert2annovar.log",
    shell:
        "{config[software][convert2annovar]} -format vcf {input} > {output} 2> {log}"


rule table_annovar:
    input:
        rules.convert2annovar.output,
    output:
        "annovar/{sample}.hg19_multianno.txt",
    benchmark:
        ".log/annovar/{sample}.table_annovar.bm"
    log:
        ".log/annovar/{sample}.table_annovar.log",
    shell:
        """
        # 获取前缀
        out_prefix=$(dirname {output})/$(basename {output} .hg19_multianno.txt)
        # 运行 table_annovar
        {config[software][table_annovar]} {input} \
        {config[database][annovar]} \
        -buildver hg19 \
        -out $out_prefix \
        -protocol refGeneWithVer,dbnsfp47a_interpro \
        -operation g,f \
        -remove -polish -nastring - \
        2> {log}
        """
