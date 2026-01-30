import pandas as pd
from pathlib import Path

diversity_files = [str(dpath) for dpath in Path(
    '/data/mengxf/Project/2025/KML250709-lvis-jiance-run1-2/result/250731/stats/').glob('*.is.diversity.tsv')]

diversity_files_remove_pos_ntc = [
    x for x in diversity_files if ('POS' not in x) and ('NTC' not in x)]

df_all = pd.DataFrame(
    columns=['sample', 'shannon', 'invsimpson', 'chao1', 'gini_coef', 'gini_simpson'])

for dfile in diversity_files_remove_pos_ntc:
    sample = Path(dfile).stem.split('.')[0]
    df = pd.read_csv(dfile, sep='\t')
    df.insert(0, 'sample', sample)
    df_all = pd.concat([df_all, df])

df_all.to_csv('diversity_summary.tsv', sep='\t', index=False)
