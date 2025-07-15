# 将汇总好的 阳控,阴控,批次控制 数据合并到原始数据中

import pandas as pd
import sys


control_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

df_control = pd.read_table(control_file, sep='\t')
df = pd.read_table(input_file, sep='\t')
df.merge(df_control, on=['Chrom', 'Start'], how='left').fillna(
    '-').to_csv(output_file, index=False, sep='\t')
