# 过滤插入位点
#   1. 1个位点 < 2 UMI (有限 UMI 数量不适用, 非随机 UMI
#   2. support reads < 10
import sys


umi_bed = sys.argv[1]
all_cov = sys.argv[2]
out_cov = sys.argv[3]

dict_stat = {}

with open(umi_bed) as f:
    for line in f:
        chrom, start, end, strand, read_count = line.strip().split('\t')[:5]
        dict_stat[(chrom, start, end)] = {'strand': strand, 'umi_num': int(read_count)}

with open(all_cov) as f:
    """
    chr10   1338137 1338201 -       1       UMI280  6       64      64      1.0000000
    chr10   42385153        42385251        +       1       UMI65   678     98      98      1.0000000
    ...
    """
    for line in f:
        lns = line.strip().split('\t')
        chrom, start, end = lns[:3]
        dict_stat[(chrom, start, end)]['all_num'] = int(lns[6])


# 输出
with open(out_cov, 'w') as f:
    f.write('#chrom\tstart\tend\tumi_num\tall_num\n')

    for k in dict_stat:
        umi_num = dict_stat[k]['umi_num']
        all_num = dict_stat[k]['all_num']

        # * 正链和负链位置不同
        chrom, start, end = k
        if dict_stat[k]['strand'] == '+':
            ostart, oend = start, str(int(start) + 1)
        else:
            ostart, oend = end, str(int(end) + 1)

        # * 过滤 support reads < 10, 过滤 UMI == 0
        if (int(all_num) < 10) or (int(umi_num) == 0):
            continue

        list_out = [chrom, ostart, oend] + [str(umi_num), str(all_num)]
        f.write('\t'.join(list_out) + '\n')
