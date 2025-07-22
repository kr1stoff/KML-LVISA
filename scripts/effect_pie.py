import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


def main(df, out_png):
    effect_count = df.groupby('Effect')['Effect'].count().sort_values(ascending=False)
    names = effect_count.index.tolist()
    # 注释到的功能变多了, 所以需要仅展示前几个
    top_num = 5
    top_counts = list(effect_count[:top_num]) + [effect_count[top_num:].sum()]
    top_names = names[:top_num] + ['Other']
    # top_counts, top_names
    plt.pie(top_counts, labels=top_names, labeldistance=1.15,
                    wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    # 画图
    ax = plt.pie(top_counts, labels=top_names, labeldistance=1.15,
                 wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    plt.savefig(out_png, dpi=300)


# * IO
combine_tsv = sys.argv[1]
out_png = sys.argv[2]

# * main
df = pd.read_table(combine_tsv, sep='\t', usecols=['Effect'])
if df.empty:
    os.system(f"touch {out_png}")
else:
    main(df, out_png)
