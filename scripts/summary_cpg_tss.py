import pandas as pd
import sys
from pathlib import Path


annotabs = sys.argv[3:]
cpg_outfile = sys.argv[1]
tss_outfile = sys.argv[2]

# annotabs = list(Path('/data/mengxf/Project/KML250703-lvis-pipeline-update/result/250714/anno').glob('*.tsv'))
annotabs = [x for x in annotabs if ('POS' not in str(x)) and ('NTC' not in str(x))]

# CpG
df_all = pd.DataFrame(columns=['Sample', 'CpG1KB', 'CpG2.5KB', 'CpG5KB', 'CpG10KB'])
for annotab in annotabs:
    df = pd.read_csv(annotab, sep='\t', usecols=['CpG1KB', 'CpG2.5KB', 'CpG5KB', 'CpG10KB'])
    df['Sample'] = Path(annotab).stem.split('.')[0]
    df_all = pd.concat([df_all, df])
# 筛选 CpG1KB 和 CpG2.5KB 列，统计非 '-' 的数量
result = (
    df_all.groupby('Sample')
    .apply(lambda x: pd.Series({
        'CpG1KB': (x['CpG1KB'] != '-').sum(),
        'CpG2.5KB': (x['CpG2.5KB'] != '-').sum(),
        'CpG5KB': (x['CpG5KB'] != '-').sum(),
        'CpG10KB': (x['CpG10KB'] != '-').sum(),
        'Unsure': len(x),
    }), include_groups=False)
)
percentage_df = result.apply(lambda x: x/x.sum(), axis=1)
percentage_df.drop(columns=['Unsure'], inplace=True)
percentage_df = percentage_df.map(lambda x: f'{x:.2%}')
percentage_df.to_csv(cpg_outfile, sep='\t')


# TSS
df_all = pd.DataFrame(columns=['Sample', 'TSS1KB', 'TSS2.5KB', 'TSS5KB', 'TSS10KB'])
for annotab in annotabs:
    df = pd.read_csv(annotab, sep='\t', usecols=['TSS1KB', 'TSS2.5KB', 'TSS5KB', 'TSS10KB'])
    df['Sample'] = Path(annotab).stem.split('.')[0]
    df_all = pd.concat([df_all, df])
# 筛选 TSS1KB 和 TSS2.5KB 列，统计非 '-' 的数量
result = (
    df_all.groupby('Sample')
    .apply(lambda x: pd.Series({
        'TSS1KB': (x['TSS1KB'] != '-').sum(),
        'TSS2.5KB': (x['TSS2.5KB'] != '-').sum(),
        'TSS5KB': (x['TSS5KB'] != '-').sum(),
        'TSS10KB': (x['TSS10KB'] != '-').sum(),
        'Unsure': len(x),
    }), include_groups=False)
)
# result.to_csv(tss_outfile, index=False, sep='\t')
percentage_df = result.apply(lambda x: x/x.sum(), axis=1)
percentage_df.drop(columns=['Unsure'], inplace=True)
percentage_df = percentage_df.map(lambda x: f'{x:.2%}')
percentage_df.to_csv(tss_outfile, sep='\t')
