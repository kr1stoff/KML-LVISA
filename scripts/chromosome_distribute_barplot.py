import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys
import os


def main(hg19_comp, df, out_png):
    # 染色体长度列表
    chrom_len_dict = {}

    with open(hg19_comp) as f:
        for line in f:
            lns = line.strip().split('\t')
            chrom_len_dict[lns[0]] = int(lns[1])

    # 根据染色体长度列表计算标准化深度
    df2 = pd.DataFrame(df.groupby('Chrom')['Depth'].sum())

    for idx in df2.index:
        df2.loc[idx, 'depth_stdz'] = df2.loc[idx, 'Depth'] / chrom_len_dict[idx]

    df2['depth_stdz_perc'] = df2['depth_stdz'] / df2['depth_stdz'].sum() * 100
    df2.reset_index(inplace=True)

    # 画图
    # x 轴顺序
    xlabel_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']

    if 'chrY' not in df2['Chrom'].tolist():
        xlabel_order.remove('chrY')

    # seaborn
    plt.figure(figsize=(15, 6))

    sns.barplot(
        df2, x="Chrom", y="depth_stdz_perc",
        native_scale=True,
        estimator="sum", errorbar=None,
        order=xlabel_order,
    )

    plt.title('Chromosome Integration Site Distribution')
    plt.xlabel('Chromosome')
    plt.ylabel('Depth Standardized Percentage (%)')
    plt.savefig(out_png, dpi=300)


# * IO
hg19_comp = sys.argv[1]  # 'script/hg19.fa.comp'
combine_tsv = sys.argv[2]   # 'anno/SRR17348516.is.combine.tsv'
out_png = sys.argv[3]   # 'stats/SRR17348516.chrom_dist.png'

# * main
df = pd.read_csv(combine_tsv, sep='\t', usecols=['Chrom', 'Depth'])
if df.empty:
    os.system(f"touch {out_png}")
else:
    main(hg19_comp, df, out_png)
