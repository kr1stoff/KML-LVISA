# 从 annovar 提取转录本号(参考或最长), 在转录本上的外显子号和 Interpro Domain 结果输出
# 按照 CHROM  POS  Transcript|Exon|Domain 格式

import sys

sys.stderr = open(snakemake.log[0], 'w')

input_file = snakemake.input[0]
output_file = snakemake.output[0]
refseq_longest_tx_file = snakemake.params[0]
