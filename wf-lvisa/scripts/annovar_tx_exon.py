# 从 annovar 提取转录本号(参考或最长), 在转录本上的外显子号
# 按照 CHROM  POS  Transcript|Exon 格式

import sys
import os
import pandas as pd

sys.stderr = open(snakemake.log[0], 'w')

input_file = snakemake.input[0]
output_file = snakemake.output[0]
refseq_longest_tx_file = snakemake.params[0]

# 如果 input_file 是空文件, 则创建空 output_file
if os.path.getsize(input_file) == 0:
    with open(output_file, 'w') as f:
        pass
else:
    # 代表 / 最长转录本字典
    df_tx = pd.read_table(refseq_longest_tx_file, sep='\t', index_col=0)
    dict_tx = df_tx.to_dict()

    # * 坐标 VCF(1-based) 转回 BED(0-based)
    # CHROM  POS  Transcript|Exon 格式
    results = []
    df = pd.read_table(input_file, sep='\t', usecols=[
        'Chr', 'Start', 'Gene.refGeneWithVer', 'AAChange.refGeneWithVer'])
    for _, row in df.iterrows():
        chrom = row['Chr']
        start = str(int(row['Start'])-1)
        gene = row['Gene.refGeneWithVer']
        aachange = row['AAChange.refGeneWithVer']
        # 转录本(参考或最长) + 外显子号
        # 如果是 intergenic 基因间功能类型, aachange 为 -
        # 有时 aachange 可能是 UNKNOWN
        if aachange == '-' or aachange == 'UNKNOWN':
            transcript, exon_num = '-', '-'
        else:
            refseq_tx = dict_tx['refseq_tx_id'].get(gene, None)
            longest_tx = dict_tx['longest_tx_id'].get(gene, None)
            aa_items = aachange.split(',')
            # 只有一个或者没找到代表和最长转录本则使用初始化的值
            transcript = aa_items[0].split(':')[1]
            exon_num = aa_items[0].split(':')[2]
            if len(aa_items) > 1:
                for aa in aa_items:
                    aas = aa.split(':')
                    aas_transcript, aa_exon_num = aas[1], aas[2]
                    if aas_transcript.split('.')[0] == longest_tx:
                        transcript, exon_num = aas_transcript, aa_exon_num
                    if aas_transcript.split('.')[0] == refseq_tx:
                        transcript, exon_num = aas_transcript, aa_exon_num
                        break
        results.append([chrom, start, f'{transcript}|{exon_num}'])
    outdf = pd.DataFrame(results)
    outdf.to_csv(output_file, sep='\t', index=False, header=False)
