# 所有样本汇总整合位点在基因的分布
# - 基因间和内含子区占比
# - 外显子区占比

import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path


# * IO
annotabs = sys.argv[4:]
effect_stat_outfile = sys.argv[1]
effect_stat_outfig = sys.argv[2]
effect_summary_file = sys.argv[3]

# * MAIN
effs = [eff for at in annotabs
        for eff in pd.read_csv(at, sep="\t", usecols=["Effect"])["Effect"].tolist()]
ser = pd.Series(Counter(effs))

# 输出统计文字
with open(effect_stat_outfile, "w") as f:
    f.write(
        f"Intron + Intergenic: {format((ser['Intron'] + ser['Intergenic']) / ser.sum(), '.2%')}\n")
    f.write(f"Exon: {format(ser['Exon'] / ser.sum(), '.2%')}\n")

# 输出统计图
ax = plt.pie(ser, labels=ser.index.tolist(), labeldistance=1.15,
             wedgeprops={"linewidth": 1, "edgecolor": "white"})
plt.savefig(effect_stat_outfig, dpi=300)

# 汇总 Effect, 每行一个样本
annotabs = [x for x in annotabs if ('POS' not in x) and ('NTC' not in x)]
df_all = pd.DataFrame(columns=['Sample', 'Effect'])
for annotab in annotabs:
    df = pd.read_csv(annotab, sep='\t', usecols=['Effect'])
    df['Sample'] = Path(annotab).stem.split('.')[0]
    df_all = pd.concat([df_all, df])
df_grpby = df_all.groupby(['Sample', 'Effect']).size().reset_index(name='Counts')
df_effect_summary = df_grpby.pivot(index='Sample', columns='Effect', values='Counts').fillna(
    0).astype(int).rename_axis(None, axis=1)
df_effect_summary.to_csv(effect_summary_file, sep='\t')
