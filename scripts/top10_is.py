# 获取 top 10 的克隆，并画图
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os


def main(df, top10_tsv, top10_png):
    # 计算 Top10
    dict_top10 = {}

    for row in df[:10].iterrows():
        loci = f"{row[1]['Chrom']}:{row[1]['Start']}_{row[1]['Gene']}"
        depth = row[1]['Depth']
        dict_top10[loci] = int(depth)

    # dict_top10['Others'] = sum(df[10:]['Depth'])
    df_top10 = pd.DataFrame.from_dict(dict_top10, orient='index', columns=['Depth'])
    df_top10.sort_values(by='Depth', inplace=True, ascending=False)
    df_others = pd.DataFrame([[sum(df[10:]['Depth'])]], columns=['Depth'], index=['Others'])
    df_top10 = pd.concat([df_top10, df_others], axis=0)
    df_top10.index.name = 'Loci'
    df_top10 = df_top10.assign(Frequency=df_top10['Depth'] / df_top10['Depth'].sum() * 100)

    # 输出表格
    df_top10.to_csv(top10_tsv, sep='\t')
    # print(df_top10)

    # 画图
    plt.figure(figsize=(10, 10))
    ticklabels = [i.split('_')[0] for i in df_top10.index]
    sns.set_style('whitegrid')
    ax = sns.barplot(data=df_top10, x=df_top10.index, y='Frequency')
    ax.set_xticklabels(ticklabels, fontdict={'rotation': 30})
    ax.set_ylabel('Frequency (%)')
    ax.set_xlabel('')
    plt.savefig(top10_png, dpi=300)


# * IO
combine_tsv = sys.argv[1]   # 'anno/SRR17348516.is.combine.tsv'
top10_tsv = sys.argv[2]  # 'stats/SRR17348516.top10_is.tsv'
top10_png = sys.argv[3]  # 'stats/SRR17348516.top10_is.png'

# * main
df = pd.read_csv(combine_tsv, sep='\t', usecols=['Chrom', 'Start', 'Depth', 'Gene'])

if df.empty:
    os.system(f"touch {top10_tsv}")
    os.system(f"touch {top10_png}")
else:
    main(df, top10_tsv, top10_png)
