"""`
去除单条非特异 read，提升整合位点识别精度。
- + 正链：左端位置频数
- - 负链：右端位置频数

设置频数阈值，去除频数小于阈值的 read (PE 双端测序需要 x2)
"""

import sys
import csv
from collections import defaultdict

# raw_bed = 'ZQX4.raw.bed'
# filtered_raw_bed = 'ZQX4.filtered.raw.bed'
# waste_raw_bed = 'ZQX4.waste.raw.bed'
# read_threas = 1
raw_bed = sys.argv[1]
filtered_raw_bed = sys.argv[2]
waste_raw_bed = sys.argv[3]
read_threas = int(sys.argv[4])

# ! 核心阈值
read_threas_pe = read_threas * 2

# 粗略的整合位点位置计数字典
raw_isite_count_dict = defaultdict(int)

# Read and count integration sites
with open(raw_bed, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for chrom, start, end, name, score, strand in reader:
        if strand == '+':
            raw_isite_count_dict[(chrom, start)] += 1
        else:
            raw_isite_count_dict[(chrom, end)] += 1

# * 过滤read数低于阈值的整合位点
remained_raw_isites = {isite for isite, count in raw_isite_count_dict.items()
                       if count > read_threas_pe}

# 重新扫描文件并分类写入, csv 库速度更快
with open(raw_bed, newline='') as f, open(filtered_raw_bed, 'w', newline='') as g, open(waste_raw_bed, 'w', newline='') as h:
    reader = csv.reader(f, delimiter='\t')
    writer_filtered = csv.writer(g, delimiter='\t')
    writer_waste = csv.writer(h, delimiter='\t')

    for line in reader:
        chrom, start, end, name, score, strand = line
        if (strand == '+' and (chrom, start) in remained_raw_isites) or (strand == '-' and (chrom, end) in remained_raw_isites):
            writer_filtered.writerow(line)
        else:
            writer_waste.writerow(line)
