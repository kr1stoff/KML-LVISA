# 所有样本汇总整合位点在基因的分布
# - 基因间和内含子区占比
# - 外显子区占比

import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


# * IO
annotabs = sys.argv[3:]
eff_stat_outfile = sys.argv[1]
eff_stat_outfig = sys.argv[2]

# * MAIN
effs = [eff for at in annotabs
        for eff in pd.read_csv(at, sep="\t", usecols=["Effect"])["Effect"].tolist()]
ser = pd.Series(Counter(effs))

# 输出统计文字
with open(eff_stat_outfile, "w") as f:
    f.write(f"Intron + Intergenic: {format((ser['Intron'] + ser['Intergenic']) / ser.sum(), '.2%')}\n")
    f.write(f"Exon: {format(ser['Exon'] / ser.sum(), '.2%')}\n")

# 输出统计图
ax = plt.pie(ser, labels=ser.index.tolist(), labeldistance=1.15,
             wedgeprops={"linewidth": 1, "edgecolor": "white"})
plt.savefig(eff_stat_outfig, dpi=300)
