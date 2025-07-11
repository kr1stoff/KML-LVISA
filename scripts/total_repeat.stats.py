# 所有样本汇总整合位点在 repeat 区域的分布, SINE/LINE/LTR/DNA

import sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


# * IO
annotabs = sys.argv[3:]
rep_stat_outfile = sys.argv[1]
rep_stat_outfig = sys.argv[2]

# * MAIN
reps = [rep for at in annotabs for rep in pd.read_csv(
    at, sep="\t", usecols=["RepClass"])["RepClass"].tolist()]
ser = pd.Series(Counter(reps))
ser.rename({"-": "UNKNOWN"}, inplace=True)
ser.sort_values(ascending=False, inplace=True)

# 不一定有 Simple_repeat
if "Simple_repeat" not in ser.index:
    ser["Simple_repeat"] = 0

# 统计
with open(rep_stat_outfile, "w") as f:
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
plt.savefig(rep_stat_outfig, dpi=300)
