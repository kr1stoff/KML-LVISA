# KML-LVISA

慢病毒整合位点分析流程

## 命令行

1. 主脚本

    ```bash
    poetry -C /data/mengxf/GitHub/KML-LVISA run python /data/mengxf/GitHub/KML-LVISA/main.py -s templates/input.tsv -w 241105
    ```

2. 单独跑 snakemake

    ```bash
    snakemake -c 32 --use-conda -s /data/mengxf/GitHub/KML-LVISA/wf-lvisa/Snakefile --configfile .temp/snakemake.yaml
    ```

## 更新

- FIXME [241226] 预计版本号 0.2.0
  - UMI 计算方法重写
  - 根据比对方向确定插入位点精确位置
  - 插入位点位置附近根据深度判断精确位置
