# 提取 annovar 结果中 interpro domains 基因功能结构域信息
# 按照 CHROM  POS  Domain 格式
import sys
import os
import pandas as pd

sys.stderr = open(snakemake.log[0], 'w')

input_file = snakemake.input[0]
output_file = snakemake.output[0]


def format_domain(domain_info):
    """格式化功能域信息，去除重复项和无效项"""
    domain_info = domain_info.replace(r'\x3b', ';').replace(r'\x2c', ',')
    items = {item for item in domain_info.split(';') if item.strip() and item.strip() != '.'}
    return ';'.join(items) if items else '-'


# 如果 input_file 是空文件, 则创建空 output_file
if os.path.getsize(input_file) == 0:
    with open(output_file, 'w') as f:
        pass
else:
    # CHROM  POS  Domain 格式
    df = pd.read_table(input_file, sep='\t', usecols=['Chr', 'Start', 'Interpro_domain'])
    # 基因结构功能域
    df['Interpro_domain'] = df['Interpro_domain'].str.replace(
        r'\x3b', ';').str.replace(r'\x2c', ',')
    df['Domain'] = df['Interpro_domain'].apply(format_domain)
    # * 坐标 VCF(1-based) 转回 BED(0-based)
    df['Start'] = df['Start'] - 1
    df[['Chr', 'Start', 'Domain']].to_csv(output_file, sep='\t', index=False, header=False)
