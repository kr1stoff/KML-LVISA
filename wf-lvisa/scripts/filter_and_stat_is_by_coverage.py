# 过滤插入位点
#   1. 1个位点 < 2 UMI (有限 UMI 数量不适用, 非随机 UMI)
#   2. support reads < 10
import sys

sys.stderr = open(snakemake.log[0], 'w')

umi_bed = snakemake.input[0]
all_cov = snakemake.input[1]
rmdup_cov = snakemake.input[2]
out_cov = snakemake.output[0]

# 统计字典
dict_stat = {}
# UMI bed. 位置信息, 链方向, UMI 数量
with open(umi_bed) as f:
    for line in f:
        chrom, start, end, strand, read_count = line.strip().split('\t')[:5]
        dict_stat[(chrom, start, end)] = {'strand': strand, 'umi_num': int(read_count)}

# 总覆盖深度文件. 位置信息, 总覆盖深度
with open(all_cov) as f:
    contents = f.readlines()
# [250731] 每个位点覆盖深度百分比
# 深度 + 深度百分比
total_depth = sum([int(line.strip().split('\t')[6]) for line in contents])
for line in contents:
    lns = line.strip().split('\t')
    chrom, start, end = lns[:3]
    dict_stat[(chrom, start, end)].update(
        {'all_num': int(lns[6]), 'all_freq': round(int(lns[6]) / total_depth, 4)})

# 过滤后的覆盖深度文件. 位置信息, 覆盖深度
with open(rmdup_cov) as f:
    for line in f:
        lns = line.strip().split('\t')
        chrom, start, end = lns[:3]
        dict_stat[(chrom, start, end)]['rmdup_num'] = int(lns[6])

# 输出
with open(out_cov, 'w') as f:
    f.write('#chrom\tstart\tend\tumi_num\tall_num\tall_freq\trmdup_num\n')
    for k in dict_stat:
        umi_num = dict_stat[k]['umi_num']
        all_num = dict_stat[k]['all_num']
        # * 正链和负链位置不同
        chrom, start, end = k
        if dict_stat[k]['strand'] == '+':
            ostart, oend = start, str(int(start) + 1)
        else:
            ostart, oend = end, str(int(end) + 1)
        # * 过滤 support reads < 10, 过滤 UMI == 0 (随机 UMI > 2)
        if (int(all_num) < 10) or (int(umi_num) == 0):
            continue
        list_out = [chrom, ostart, oend, str(umi_num), str(all_num),
                    str(dict_stat[k]['all_freq']), str(dict_stat[k]['rmdup_num'])]
        f.write('\t'.join(list_out) + '\n')
