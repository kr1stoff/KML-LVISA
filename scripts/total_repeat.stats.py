# 所有样本汇总整合位点在 repeat 区域的分布, SINE/LINE/LTR/DNA

import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path


# * IO
annotabs = sys.argv[4:]
repeat_stat_outfile = sys.argv[1]
repeat_stat_outfig = sys.argv[2]
repclass_summmary_outfile = sys.argv[3]

# * MAIN
reps = [rep for at in annotabs
        for rep in pd.read_csv(at, sep="\t", usecols=["RepClass"])["RepClass"].tolist()]
ser = pd.Series(Counter(reps))
ser.rename({"-": "UNKNOWN"}, inplace=True)
ser.sort_values(ascending=False, inplace=True)

# 不一定有 Simple_repeat
if "Simple_repeat" not in ser.index:
    ser["Simple_repeat"] = 0

# 统计
with open(repeat_stat_outfile, "w") as f:
    f.write(f'SINE: {format((ser["SINE"]) / ser.sum(), ".2%")}\n')
    f.write(f'LINE: {format(ser["LINE"] / ser.sum(), ".2%")}\n')
    f.write(f'LTR: {format(ser["LTR"] / ser.sum(), ".2%")}\n')
    f.write(f'DNA: {format(ser["DNA"] / ser.sum(), ".2%")}\n')
    f.write(f'Simple_repeat: {format(ser["Simple_repeat"] / ser.sum(), ".2%")}\n')

# 画图
names = ser.index.tolist()
newnms = [names[i] if i < 5 else "" for i in range(len(names))]
ax = plt.pie(ser, labels=newnms, labeldistance=1.15,
             wedgeprops={"linewidth": 1, "edgecolor": "white"})
plt.savefig(repeat_stat_outfig, dpi=300)

# 汇总
annotabs = [x for x in annotabs if ('POS' not in x) and ('NTC' not in x)]
df_all = pd.DataFrame(columns=['Sample', 'RepClass'])
for annotab in annotabs:
    df = pd.read_csv(annotab, sep='\t', usecols=['RepClass'])
    df['Sample'] = Path(annotab).stem.split('.')[0]
    df_all = pd.concat([df_all, df])
    df_all = df_all[df_all['RepClass'] != '-']
df_grpby = df_all.groupby(['Sample', 'RepClass']).size().reset_index(name='Counts')
df_repclass_summary = df_grpby.pivot(index='Sample', columns='RepClass', values='Counts').fillna(
    0).astype(int).rename_axis(None, axis=1)
df_repclass_summary.to_csv(repclass_summmary_outfile, sep='\t')
