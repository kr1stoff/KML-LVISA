# KML-LVISA

慢病毒整合位点分析流程

## 命令行

1. 主脚本

    ```bash
    # 安装了 poetry 的环境
    conda activate python3.10
    poetry -C /data/mengxf/GitHub/KML-LVISA run python /data/mengxf/GitHub/KML-LVISA/main.py -s templates/input.tsv -w 241105
    ```

2. 单独跑 snakemake

    ```bash
    # 安装了 snakemake 的环境
    snakemake -c 32 --use-conda -s /data/mengxf/GitHub/KML-LVISA/wf-lvisa/Snakefile --configfile .temp/snakemake.yaml
    ```

## 更新

- [241226] 版本号 0.2.0
  - UMI 计算方法重写, 加入艾吉泰康 384 UMI 文件
  - 根据比对方向确定插入位点精确位置
