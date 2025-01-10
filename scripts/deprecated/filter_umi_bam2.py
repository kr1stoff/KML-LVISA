import sys
import json
import Levenshtein
import re


def parse_sam_line(line) -> tuple:
    """
    解析 SAM 文件行
    :param line: SAM 文件行
    :return: (umi_id, rname, start, end, strand, tlen, read_len)
    """
    lns = line.strip().split('\t')
    qname, flag, rname, pos, _, cigar = lns[:6]
    tlen = abs(int(lns[8]))
    read_len = int(re.search('(\d+)M', cigar).group(1))
    umi_seq = qname.strip().split(':')[-1]
    strand = '-' if int(flag) & 16 else '+'
    # 位置
    start = int(pos)
    end = int(start) + read_len
    # 负链, 交换位置
    if strand == '-':
        start, end = end, start
    umi_id = assign_umi_seq_to_id(umi_seq)
    return (umi_id, rname, start, end, strand, tlen, read_len, flag)


def assign_umi_seq_to_id(umi_seq):
    """根据 UMI 序列计算 UMI ID"""
    for k, v in umi_dict.items():
        if Levenshtein.distance(umi_seq, k) <= 1:
            return v
    return 'UMI0'


def add_umi_to_list(bed_dict, counter, umi_id):
    """重复步骤，添加 UMI 到列表"""
    bed_dict[counter].setdefault('umi', [])
    if umi_id not in bed_dict[counter]['umi']:
        bed_dict[counter]['umi'].append(umi_id)


def sam_to_bed_dict() -> dict:
    # SAM
    # * Fields
    # 1. QNAME read 名, UMI 序列
    # 2. FLAG 看方向
    # 3. RNAME 参考序列名
    # 4. POS 位置
    # 6. CIGAR 比对信息, 用 M 和 序列方向给定起始和终止位置
    # 9. TLEN 两个 reads 之间的距离，根据距离判断文库大小，使用阈值 1000 过滤
    # * 筛选规则
    # 1. 比对的长度大于 50bp 虽然前面 samtools 已经过滤了，但是还是再过滤一次
    # 2. 两个 reads 之间的距离小于 1000

    bed_dict = {}
    counter = 0
    with open(in_sam) as f:

        while True:
            """确定初始行"""
            line = next(f)
            umi_id, rname_1, start_1, end_1, strand, tlen, read_len, flag = parse_sam_line(line)
            # 跳过 read2
            if int(flag) & 128:
                continue
            # 满足条件则终止循环
            if (read_len > 50) and (tlen < 1000):
                break
        # 初始字典信息，包括位置、方向、UMI
        bed_dict[counter] = {'postion': (rname_1, start_1, end_1),
                             'strand': strand, 'umi': [umi_id]}

        for line in f:
            umi_id, rname, start, end, strand, tlen, read_len, flag = parse_sam_line(line)
            # 过滤 + 跳过 read2
            if (read_len < 50) or (tlen > 1000) or (int(flag) & 128):
                continue

            # 如果染色提相同
            if rname == rname_1:
                # 如果新的一行包含在之前的区域
                if (start <= start_1) and (end >= end_1):
                    bed_dict[counter]['postion'] = (rname, start, end)
                    add_umi_to_list(bed_dict, counter, umi_id)
                # 如果新的一行与之前区域 overlap 左
                elif (start <= start_1) and (start_1 <= end <= end_1):
                    bed_dict[counter]['postion'] = (rname, start, end_1)
                    add_umi_to_list(bed_dict, counter, umi_id)
                # 如果新的一行与之前区域 overlap 右
                elif (end_1 >= start >= start_1) and (end >= end_1):
                    bed_dict[counter]['postion'] = (rname, start_1, end)
                    add_umi_to_list(bed_dict, counter, umi_id)
                else:
                    counter += 1
                    bed_dict[counter] = {'postion': (rname, start, end),
                                         'strand': strand, 'umi': [umi_id]}
                    rname_1, start_1, end_1 = rname, start, end
            else:
                counter += 1
                bed_dict[counter] = {'postion': (rname, start, end),
                                     'strand': strand, 'umi': [umi_id]}
                rname_1, start_1, end_1 = rname, start, end

        return bed_dict


def output_bed_file(bed_dict):
    # 输出bed文件
    with open(out_bed, 'w') as f:
        for i in bed_dict:
            position = list(map(str, bed_dict[i]['postion']))
            strand = bed_dict[i]['strand']
            if strand == '-':
                position[1], position[2] = position[2], position[1]
            umi = ','.join(bed_dict[i]['umi'])
            umi_num = str(len(list(filter(lambda x: x != 'UMI0', bed_dict[i]['umi']))))
            f.write('\t'.join(position + [strand, umi_num, umi]) + '\n')


# main
# * 需要 samtools sort 过的 SAM 文件
in_sam = sys.argv[1]
ajtk_umi = sys.argv[2]
out_bed = sys.argv[3]
# * UMI字典, 经验证 UMI 间编辑距离最小为 3
umi_dict = json.load(open(ajtk_umi, 'r'))

bed_dict = sam_to_bed_dict()
output_bed_file(bed_dict)
