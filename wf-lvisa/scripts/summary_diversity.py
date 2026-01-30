import pandas as pd
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], "w")

infiles = snakemake.input
outfile = snakemake.output[0]

infiles_remove_pos_ntc = [x for x in infiles
                          if ('POS' not in x) and ('NTC' not in x)]

df_all = pd.DataFrame(
    columns=['sample', 'shannon', 'invsimpson', 'chao1', 'gini_coef', 'gini_simpson'])

for dfile in infiles_remove_pos_ntc:
    sample = Path(dfile).stem.split('.')[0]
    df = pd.read_csv(dfile, sep='\t')
    df.insert(0, 'sample', sample)
    df_all = pd.concat([df_all, df])

df_all.to_csv(outfile, sep='\t', index=False)
