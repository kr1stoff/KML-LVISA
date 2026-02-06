# KML-LVISA

![GitHub followers](https://img.shields.io/github/followers/kr1stoff)
![GitHub Created At](https://img.shields.io/github/created-at/kr1stoff/KML-LVISA)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/kr1stoff/KML-LVISA)
![GitHub Release](https://img.shields.io/github/v/release/kr1stoff/KML-LVISA)
![GitHub License](https://img.shields.io/github/license/kr1stoff/KML-LVISA)

慢病毒整合位点分析流程

## 命令行

1. 拆分数据

    ```bash
    bcl2fastq --no-lane-splitting --barcode-mismatches 1 --processing-threads 32 \
        --runfolder-dir /data/rawdata/illumina/NEXTseq500/250707_NB501947_0941_AHKKYYBGXW \
        --input-dir /data/rawdata/illumina/NEXTseq500/250707_NB501947_0941_AHKKYYBGXW/Data/Intensities/BaseCalls \
        --sample-sheet samplesheet.csv \
        --output-dir FASTQ
    ```

2. 准备流程输入文件

    ```bash
    ~/miniforge3/envs/python3.12/bin/python /data/mengxf/GitHub/KML-LVISA/utils/prepare_input_from_bcl2fastq_out.py \
        --input-dir /data/mengxf/Project/KML250620-lvis-jiance-run1/FASTQ \
        --output-file /data/mengxf/Project/KML250620-lvis-jiance-run1/input.tsv
    ```

3. 主脚本
  使用绝对路径

    ```bash
    ~/miniforge3/envs/python3.12/bin/python /data/mengxf/GitHub/KML-LVISA/main.py \
      --sample-table /data/mengxf/Project/KML250513_lvisa_update/input.tsv \
      --work-dir /data/mengxf/Project/KML250513_lvisa_update/result/250514
    ```

4. (可选)单独跑 snakemake

    ```bash
    # 安装了 snakemake 的环境
    snakemake -c 32 --use-conda -s /data/mengxf/GitHub/KML-LVISA/wf-lvisa/Snakefile --configfile .temp/snakemake.yaml --scheduler greedy --ignore-incomplete
    ```

5. 同步目录

    不同步大文件 FASTQ, BAM 等

    ```bash
    rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' --exclude '3ltr/' \
      /data/mengxf/Project/KML250709-lvis-jiance-run1-2/result/250731/ \
      /data/share/samba/public/bioinfo/KML250709-lvis-jiance-run1-2/result/250731/
    ```

    同步 SAV 文件

    ```bash
    rsync -auvP --delete \
      --include 'Images/***' --include 'InterOp/***' --include 'Thumbnail_Images/***' \
      --include 'RunInfo.xml' --include 'RunParameters.xml' \
      --exclude '*' \
      /data/rawdata/illumina/NEXTseq500/250707_NB501947_0941_AHKKYYBGXW/ \
      /data/share/samba/public/bioinfo/KML250709-lvis-jiance-run1-2/250707_NB501947_0941_AHKKYYBGXW/
    ```

## 说明

## 更新

- [20260204] 0.5.8
  - fastp 切不干净, 加一步 cutadapt. UMI前的一段接头R1用正向序列 ^TGGATAAAGTCGGA 切3'端', R2 用反向互补序列 TCCGACTTTATCCA$ 切5'端'
  - 新增 57bp LTR 序列判读, 文库上内引物到整合位点全部的 57bp LTR 序列, 实现更精准判读 read 是否包含 LTR, 是否为有效整合位点数据

- [20260203] 0.5.7
  - fastp 调整参数, 删除 --adapter_sequence_r2 TCCGACTTTATCCA --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction; 在 fastp 不做 umi_adapter 切除, 后续 cutadapt 来切; 新增参数 --trim_front2 14, 切 R2 UMI 拆分后的5'端14bp, 默认这部分是 umi_adapater

- [20260130] 0.5.6
  - 新增 diversity 汇总文件

- [20251223] 0.5.5
  - 修复 BUG. 先过滤后去重, reads 不成对, 没有 ms(mate score) tag, samtools markdup 会报错. 调整为先去重后过滤(该步骤与仅影响去重reads计数).
  - 新增 clean reads/bases 占比
  - 优化 map2hg19 部分日志, 管道符每条命令都要 2>>{log}

- [250815] 0.5.4
  - 新增过滤整合位点脚本
    - Batch <= 0.9
    - UMIs >= 3
    - 10 <= Depth/UMI <= 100 (384 全保留)
    - 删除 RmdupDepth 列
  - TSS 换参考, 使用 UCSC RefSeq FANTOM TSS
  - CpG, TSS, effect, repeat 样本结果汇总表都用小数占比数值

- [250806] 0.5.3
  - (LYQ) 新增染色体汇总文件
  - 修复 OncoKB 汇总部分 BUG, 由于注释软件从 snpEff 转至 annovar, 外显子从 Exon 变为 exonic, 引起错误
  - 优化单样本染色体分布图

- [250730] 0.5.2
  - 单样本注释结果中新增
    - 插入位点 reads 占比 (参考 TCR)
    - 插入位点去重 reads 数

- [260716] 0.5.1
  - 注释软件, snpEff -> annovar
  - (MLL) 基因功能结构域, interpro_domain
  - (MLL) 转录本号(代表性或最长) + 外显子数

- [250630] 0.5.0
  - 新增输出 fastp.stats.xlsx, 便于查看
  - (LYQ) 单样本 .is.combine.tsv 文件更新注释信息
    - cpg 新增 1kb, 2.5kb, 5kb, 10kb 距离分布
    - TSS 新增 1kb, 2.5kb, 5kb, 10kb 距离分布
    - 每 1MB 基因组区域 GC 值
    - 新增 Depth / UMI 比值列
  - (MXFA) 单样本 .is.combine.tsv 文件更新阴控阳控批次信息
    - 批次内频率, 解读时使用 50% 频率筛选
    - 阴控有无
    - 阳控有无 (与阴控类似, 所有检出的位点都列出)
  - (MXFA) BAM 过滤参数 TLEN 调整 (1000 -> 600), 要求双端 read 距离更近
  - (MXFA) UMI support reads 1 -> 3, 结果太多提高阈值
  - (MXFA) LTR参考文件 16bp 改回 25bp
  - 更新综合报告部分内容, 保留 *0.4.0* 更新内容. 汇总各样本数值到一个表格, 每行一个样本
    - effect
    - repeat
    - CpG
    - TSS

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
