# 整合位点安全性分析
# - TSGs / Oncogene
# 位点数量：外显子上的注释到的数量
# 位点信息：基因名称？坐标位置？

import sys
from pathlib import Path
import pandas as pd


# * IO
annotabs = sys.argv[2:]
onco_stat_outfile = sys.argv[1]


# * MAIN
with open(onco_stat_outfile, "w") as f:
    f.write("样本实验号\tTSGs 位点数量\tTSGs 位点信息\tOncogene 位点数量\tOncogene 位点信息\n")
    for at in annotabs:
        df = pd.read_csv(at, sep="\t", usecols=["Chrom", "Start", "Gene", "Effect", "Oncogene/TSG"])
        # 在外显子上注释到 Oncogene, 有重复时去重
        oncodf = df[(df["Effect"] == "exonic") & (df["Oncogene/TSG"] == "oncogene")]
        onco_num = oncodf.shape[0]
        onco_genes = ",".join(set(oncodf["Gene"].tolist())) if onco_num else "-"
        # 在外显子上注释到 TSG
        tsgdf = df[(df["Effect"] == "exonic") & (df["Oncogene/TSG"] == "tsg")]
        tsg_num = tsgdf.shape[0]
        tsg_genes = ",".join(set(tsgdf["Gene"].tolist())) if tsg_num else "-"
        sample = Path(at).stem.split(".")[0]
        f.write(f"{sample}\t{tsg_num}\t{tsg_genes}\t{onco_num}\t{onco_genes}\n")
