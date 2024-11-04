# 过滤插入位点
#   1. 1个位点 < 2 UMI (羿圣试剂盒不适用)
#   2. support reads < 10
import sys
import pdb


umi_cov = sys.argv[1]
all_cov = sys.argv[2]
out_cov = sys.argv[3]

dict_stat = {}

with open(umi_cov) as f:
    for line in f:
        chrom, start, end, reads = line.strip().split('\t')[:4]
        dict_stat.setdefault((chrom, start, end), {})
        dict_stat[(chrom, start, end)]['umi_num'] = int(reads)


with open(all_cov) as f:
    for line in f:
        chrom, start, end, reads = line.strip().split('\t')[:4]
        dict_stat[(chrom, start, end)]['all_num'] = int(reads)

# 暂时不排序，最后整合出再排序
# sorted_dict_stat = dict(sorted(dict_stat.items(), key=lambda it: it[1]['umi_num'], reverse=True))

# 输出
with open(out_cov, 'w') as f:
    f.write('#chrom\tstart\tend\tumi_num\tall_num\n')

    for k in dict_stat:
        umi_num = dict_stat[k]['umi_num']
        all_num = dict_stat[k]['all_num']

        # 过滤
        # loci umi < 2 (羿圣试剂盒不适用)
        # support reads < 10
        # if any([int(umi_num) < 2, int(all_num) < 10]):
        if int(all_num) < 10:
            continue

        list_out = list(k) + [str(umi_num), str(all_num)]
        f.write('\t'.join(list_out) + '\n')
