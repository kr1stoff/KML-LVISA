import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys
import os
from collections import Counter

sys.stderr = open(snakemake.log[0], 'w')


def main(hg19_comp, chroms, out_png):
    # 染色体长度
    chrom_len_dict = {}
    with open(hg19_comp) as f:
        for line in f:
            lns = line.strip().split('\t')
            if lns[0] == 'chrY':
                continue
            chrom_len_dict[lns[0]] = int(lns[1])
    # 当前样本染色体分布
    chroms = pd.read_csv(combine_tsv, sep='\t', usecols=['Chrom'])['Chrom'].tolist()
    chrom_dict = dict(Counter(chroms))
    df = pd.DataFrame([chrom_len_dict, chrom_dict]).T
    df.columns = ['length', 'count']
    df['Chromosome length'] = df['length'] / df['length'].sum() * 100
    df['Intergration sites'] = df['count'] / df['count'].sum() * 100
    df.reset_index(inplace=True, names='Chromosome')
    df_melt = df.melt(id_vars=['Chromosome'], value_vars=[
                      'Chromosome length', 'Intergration sites'])
    # 作图
    xlabel_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrM']
    plt.figure(figsize=(10, 10))
    sns.barplot(df_melt, y="Chromosome", x="value", hue="variable", order=xlabel_order)
    plt.xlabel('Intergration sites (%)')
    plt.savefig(out_png, dpi=300)


# * IO
hg19_comp = snakemake.input[1]
combine_tsv = snakemake.input[0]
out_png = snakemake.output[0]

# * main
chroms = pd.read_csv(combine_tsv, sep='\t', usecols=['Chrom'])['Chrom'].tolist()
if not chroms:
    os.system(f"touch {out_png}")
else:
    main(hg19_comp, chroms, out_png)
