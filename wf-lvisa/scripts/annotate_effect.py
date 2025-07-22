# 使用 annovar 信息注释插入位点. 从 snpEff 改用 annovar
# chromosome  start  effect|gene

import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], 'w')

input_file = snakemake.input[0]
output_file = snakemake.output[0]

# 如果 input_file 是空文件, 则创建空 output_file
if os.path.getsize(input_file) == 0:
    with open(output_file, 'w') as f:
        pass
else:
    df = pd.read_table(input_file, sep='\t', usecols=[
        'Chr', 'Start', 'Func.refGeneWithVer', 'Gene.refGeneWithVer'])
    df['effect_and_gene'] = df['Func.refGeneWithVer'] + '|' + df['Gene.refGeneWithVer']
    # * 坐标 VCF(1-based) 转回 BED(0-based)
    df['Start'] = df['Start'] - 1
    df[['Chr', 'Start', 'effect_and_gene']].to_csv(output_file, sep='\t', index=False, header=False)
