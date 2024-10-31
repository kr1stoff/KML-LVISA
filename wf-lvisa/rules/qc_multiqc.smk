rule multiqc:
    input:
        expand('qc/fastqc/{sample}', sample=config['samples'])
    output:
        directory('qc/multiqc')
    benchmark:
        '.log/multiqc/multiqc.bm'
    log:
        '.log/multiqc/multiqc.log'
    conda:
        config['conda']['basic']
    shell:
        """
        #生成交互式的可视化网页
        multiqc {input} --outdir {output} > {log} 2>&1
        """
