import os
import pandas as pd
import matplotlib.pyplot as plt
import sys


def main(df, out_png):
    repclass_count = df.groupby('Rep Class')['Rep Class'].count()
    # 合并 - unnknown other 三种注释重复类型
    # 处理 BUG - KeyError: 'Unknown'
    if 'Other' not in repclass_count:
        repclass_count['Other'] = 0
    if '-' in repclass_count:
        repclass_count['Other'] += repclass_count['-']
        repclass_count.drop(['-'], inplace=True)
    if 'Unknown' in repclass_count:
        repclass_count['Other'] += repclass_count['Unknown']
        repclass_count.drop(['Unknown'], inplace=True)

    repclass_count.sort_values(ascending=False, inplace=True)
    names = repclass_count.index.tolist()
    # 只显示前 6 个
    if len(names) > 6:
        names[6:] = [''] * (len(names) - 6)

    # 画图
    ax = plt.pie(repclass_count, labels=names, labeldistance=1.15,
                 wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    plt.savefig(out_png, dpi=300)


# * IO
combine_tsv = sys.argv[1]   # 'anno/SRR17348516.is.combine.tsv'
out_png = sys.argv[2]   # 'stats/SRR17348516.repclass_pie.png'

# * main
df = pd.read_table(combine_tsv, sep='\t', usecols=['Rep Class'])
if df.empty:
    os.system(f"touch {out_png}")
else:
    main(df, out_png)
