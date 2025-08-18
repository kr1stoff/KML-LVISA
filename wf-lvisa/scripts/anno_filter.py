# [250815 LYQ] 筛选条件
# 1. Batch <= 0.9
# 2. UMI >= 3
# 3. 10 <= Depth/UMI <= 100 (384 除外)
# 删除 RmdupDepth 列

import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], 'w')

infile = snakemake.input[0]
outfile = snakemake.output[0]

# 阴控阳控不要做筛选
if ('POS' in infile) or ('NTC' in infile):
    df = pd.read_csv(infile, sep='\t')
    dfout = df.drop(columns='RmdupDepth')
    dfout.to_csv(outfile, sep='\t', index=False)
else:
    df = pd.read_csv(infile, sep='\t')
    df1 = df[df['Batch'] <= 0.9]
    df2 = df1[df1['UMIs'] >= 3]
    df3_384 = df2[df2['UMIs'] == 384]
    df3_lt384 = df2[df2['UMIs'] < 384]
    df3_2 = df3_lt384[(df3_lt384['Depth/UMI'] >= 10.0) & (df3_lt384['Depth/UMI'] <= 100.0)]
    df3 = pd.concat([df3_384, df3_2])
    dfout = df3.drop(columns='RmdupDepth')

    dfout.to_csv(outfile, sep='\t', index=False)
