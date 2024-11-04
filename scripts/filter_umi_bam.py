# 过滤 1 个 UMI < 3 条 reads 支持
import sys
import pysam


in_bam = sys.argv[1]
out_bam = sys.argv[2]

# gencore sam tag FR:i:1  RR:i:5, FF 表示正向支持 reads 数, RR 表示反向支持 reads 数
bam = pysam.AlignmentFile(in_bam, 'rb')

with pysam.AlignmentFile(out_bam, 'wb', template=bam) as f:
    for read in bam:
        dict_tags = {t[0]: t[1] for t in read.get_tags()}
        # 正向 reads
        umi_count = int(dict_tags['FR'])
        # 反向 reads
        if 'RR' in dict_tags:
            umi_count += int(dict_tags['RR'])

        # 1 个 UMI < 3 条 reads 支持, 过滤
        if umi_count < 3:
            continue

        # 写结果 bam
        f.write(read)

bam.close()
