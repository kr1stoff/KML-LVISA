import pandas as pd
import matplotlib.pyplot as plt


combine_tsv = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.combine.tsv'
out_png = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/stats/SRR17348516.repclass_pie.png'

df = pd.read_table(combine_tsv, sep='\t', usecols=['Rep Class'])
repclass_count = df.groupby('Rep Class')['Rep Class'].count()
# 合并 - unnknown other 三种注释重复类型
repclass_count['Other'] = repclass_count['-'] + repclass_count['Unknown'] + repclass_count['Other']
repclass_count.drop(['-', 'Unknown'], inplace=True)
repclass_count.sort_values(ascending=False, inplace=True)
names = repclass_count.index.tolist()
# 只显示前 6 个
if len(names) > 6:
    names[6:] = [''] * (len(names) - 6)

# 画图
ax = plt.pie(repclass_count, labels=names, labeldistance=1.15, wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
plt.savefig(out_png, dpi=300)
