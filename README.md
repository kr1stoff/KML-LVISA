# KML-LVISA

慢病毒整合位点分析流程

## 命令行

1. 主脚本

    ```bash
    poetry -C /data/mengxf/GitHub/KML-LVISA run python main.py -s templates/input.tsv -w 241105
    ```

2. 单独跑 snakemake

    ```bash
    snakemake -c 32 --use-conda -s wf-lvisa/Snakefile --configfile .temp/snakemake.yaml
    ```
