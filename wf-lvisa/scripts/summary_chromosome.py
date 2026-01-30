import sys
from pathlib import Path
import pandas as pd


sys.stderr = open(snakemake.log[0], "w")


annotabs = snakemake.input
chrom_summmary_outfile = snakemake.output[0]

annotabs = [x for x in annotabs if ('POS' not in x) and ('NTC' not in x)]

df_all = pd.DataFrame(columns=['Sample', 'Chrom'])
for annotab in annotabs:
    df = pd.read_csv(annotab, sep='\t', usecols=['Chrom'])
    df['Sample'] = Path(annotab).stem.split('.')[0]
    df_all = pd.concat([df_all, df])
df_grpby = df_all.groupby(['Sample', 'Chrom']).size().reset_index(name='Counts')
df_grpby = df_grpby[df_grpby['Chrom'] != 'chrY']
chrom_order = [f'chr{i}' for i in range(1,23)] + ['chrX', 'chrM']
df_grpby['Chrom'] = pd.Categorical(df_grpby['Chrom'], categories=chrom_order, ordered=True)
df_grpby.sort_values('Chrom', inplace=True)
outdf = df_grpby.pivot(index='Sample', columns='Chrom', values='Counts').fillna(0).astype(
    int).rename_axis(None, axis=1)

outdf.to_csv(chrom_summmary_outfile, sep='\t')
