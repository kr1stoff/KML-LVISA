# KML-LVISA

慢病毒整合位点分析流程

## 命令行

0. 准备流程输入文件

    ```bash
    poetry run -C /data/mengxf/GitHub/KML-LVISA \
        python /data/mengxf/GitHub/KML-LVISA/scripts/prepare_input_from_bcl2fastq_out.py \
        --input-dir /data/mengxf/Project/KML250620-lvis-jiance-run1/FASTQ \
        --output-file /data/mengxf/Project/KML250620-lvis-jiance-run1/input.tsv
    ```

1. 主脚本
  使用绝对路径

    ```bash
    poetry -C /data/mengxf/GitHub/KML-LVISA run \
      python /data/mengxf/GitHub/KML-LVISA/main.py \
      --sample-table /data/mengxf/Project/KML250513_lvisa_update/input.tsv \
      --work-dir /data/mengxf/Project/KML250513_lvisa_update/result/250514
    ```

2. 单独跑 snakemake

    ```bash
    # 安装了 snakemake 的环境
    mamba run -n snakemake \
      snakemake -c 32 --use-conda \
      -s /data/mengxf/GitHub/KML-LVISA/wf-lvisa/Snakefile \
      --configfile .temp/snakemake.yaml \
      --scheduler greedy
    ```

3. 同步目录
  不同步大文件 FASTQ, BAM 等

    ```bash
    rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' --exclude '3ltr/' \
      /data/mengxf/Project/KML250513_lvisa_update/ \
      /data/share/samba/public/bioinformatics/KML250513_lvisa_update/
    ```

## 更新

- [250630] 0.4.2
  - 更新综合报告部分内容, 保留 *0.4.0* 更新内容. 汇总各样本数值到一个表格, 每行一个样本
    - effect
    - repeat
    - cpg, 新增 1kb, 2.5kb, 5kb, 10kb 距离分布
    - TSS, 新增 1kb, 2.5kb, 5kb, 10kb 距离分布
    - 插入位点在基因组不同 GC 含量区域的分布情况, 先看 UCSC 是否有数据库

- [250623] 0.4.1
  - 新增 prepare_input_from_bcl2fastq_out.py 脚本
    - 从 bcl2fastq 输出的 FASTQ 文件夹生成输入文件
    - 生成的输入文件格式为 sampleid, read1, read2 三列制表符分隔
  - 修复 snakemake 调度 BUG
    - Snakemake 使用 `--scheduler greedy`

- [250513] 0.4.0
  - 针对综合报告部分内容更新流程
    - 总体整合位点在基因分布
    - 总体整合位点在 repeat 区域分布
    - 整合位点安全性分析，TSGs/Oncogene 位点数量和信息
  - 多线程复制 FASTQ

- [250206] 0.3.1
  - 针对性能验证样本 ZQX4 单条非特异 read 噪音过滤
  - (AGTK 建议) 3LTR 引物比对位置过滤. awk '$11>30', 引物长度为 25，允许 5bp 偏移
  - UMI support reads 3 -> 1

- [250123] 0.3.0
  - 修复 NTC 无数据流程中断的 BUG

- [250120] 0.2.1
  - 添加 HGNC 基因全名. FRYL: FRY like transcription coactivator

- [241226] 0.2.0
  - UMI 计算方法重写, 加入 AGTK 384 UMI 文件
  - 根据比对方向确定插入位点精确位置

## 报错

- pulp.apis.core.PulpSolverError: Pulp: Error while trying to execute, use msg=True for more detailscbc
    修改 Snakemake 调度设置
    尝试修改 Snakemake 的调度策略，避免使用 ILP（整数线性规划）调度：

    ```bash
    snakemake --scheduler greedy
    # 或者
    snakemake --scheduler ilp
    ```
