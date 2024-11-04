import pandas as pd
import matplotlib.pyplot as plt


combine_tsv = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.combine.tsv'
out_png = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/stats/SRR17348516.effect_pie.png'

df = pd.read_table(combine_tsv, sep='\t', usecols=['Effect'])
effect_count = df.groupby('Effect')['Effect'].count()
names = effect_count.index.tolist()

# 画图
ax = plt.pie(effect_count, labels=names, labeldistance=1.15, wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
plt.savefig(out_png, dpi=300)
