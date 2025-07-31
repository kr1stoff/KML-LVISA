# 根据 anno/*.is.combine.tsv 文件，生成 control.tsv 文件.
# 阴控标识 PC, 阳控标识 NTC
# 计算批次频率排除 PC 和 NTC

import sys
import pandas as pd
from collections import defaultdict
from functools import reduce, partial


input_anno_files = sys.argv[2:]
output_control_file = sys.argv[1]


# 样本位点字典
sample_site_dict = defaultdict(int)
sample_num = 0

# 初始化阴控/阳控 DataFrame
df_ntc = pd.DataFrame(columns=['Chrom', 'Start', 'NTC'])
df_pc = pd.DataFrame(columns=['Chrom', 'Start', 'PC'])


# 迭代输入文件, 分别处理阳控、阴控和样本
for infile in input_anno_files:
    if 'POS' in infile:
        # 阳控
        df_pc = pd.read_table(infile, sep='\t', usecols=['Chrom', 'Start'], dtype=object)
        df_pc['PC'] = 'PC'
    elif 'NTC' in infile:
        # 阴控
        df_ntc = pd.read_table(infile, sep='\t', usecols=['Chrom', 'Start'], dtype=object)
        df_ntc['NTC'] = 'NTC'
    else:
        # 检测样本
        sample_num += 1
        df_batch = pd.read_table(infile, sep='\t', usecols=['Chrom', 'Start'], dtype=object)
        for _, row in df_batch.iterrows():
            sample_site_dict[(row['Chrom'], row['Start'])] += 1

# 批次频率
df_batch = pd.DataFrame([(k[0], k[1], v) for k, v in sample_site_dict.items()],
                        columns=['Chrom', 'Start', 'Count'])
df_batch['Batch'] = round(df_batch['Count'] / sample_num, 4)
df_batch.drop('Count', axis=1, inplace=True)

# 合并阳控,阴控和样本
mymerge = partial(pd.merge, on=['Chrom', 'Start'], how='outer')
df = reduce(mymerge, [df_pc, df_ntc, df_batch])
df = df.fillna('-')
df.to_csv(output_control_file, sep='\t', index=False)
