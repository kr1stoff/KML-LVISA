# KML-LVISA

慢病毒整合位点分析流程

## 命令行

1. 主脚本

    ```bash
    poetry -C /data/mengxf/GitHub/KML-LVISA run python /data/mengxf/GitHub/KML-LVISA/main.py -s templates/input.tsv -w 241105
    ```

2. 单独跑 snakemake

    ```bash
    # 安装了 snakemake 的环境
    snakemake -c 32 --use-conda -s /data/mengxf/GitHub/KML-LVISA/wf-lvisa/Snakefile --configfile .temp/snakemake.yaml
    ```

3. 同步目录
    - 不同步大文件 FASTQ, BAM 等

        ```bash
        rsync -auvP --delete --include 'anno/***' --include 'stats/***' --exclude '*fastq.gz' --include 'qc/***' --exclude '*' /data/mengxf/Project/KML250214_lvis_YanZhengRUN3/results/250214/ /data/share/samba/public/bioinformatics/250212_NB501947_0938_AHH2WMAFX7/250214/
        ```

## 更新

- [250206] 0.3.1
  - 针对性能验证样本 ZQX4 单条非特异 read 噪音过滤
  - (AGTK建议) 3LTR 引物比对位置过滤. awk '$11>30', 引物长度为 25，允许 5bp 偏移
  - UMI support reads 3 -> 1

- [250123] 0.3.0
  - 修复 NTC 无数据流程中断的 BUG

- [250120] 0.2.1
  - 添加 HGNC 基因全名. FRYL: FRY like transcription coactivator

- [241226] 0.2.0
  - UMI 计算方法重写, 加入 AGTK 384 UMI 文件
  - 根据比对方向确定插入位点精确位置
