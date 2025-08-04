# 合并统计结果到一个表格中

import sys


# 输入文件参数
is_cov_file, effect_file, oncokb_file, \
    cpg1kb_file, cpg2d5kb_file, cpg5kb_file, cpg10kb_file, \
    tss1kb_file, tss2d5kb_file, tss5kb_file, tss10kb_file, \
    rmsk_file, gc1mb_file, fullname_file, \
    tx_exon_file, domain_file, \
    combine_out = sys.argv[1:]


def get_pos_anno_dict(anno_file):
    """
    获取注释信息.
    :param anno_file: 注释文件
    :return: 注释字典 {(chrom, start): anno_info}
    """
    curr_anno_dict = {}
    with open(anno_file) as f:
        for line in f:
            chrom, start, anno_info = line.strip().split('\t')[:3]
            # * 重复的位置只要第一个
            if (chrom, start) not in curr_anno_dict:
                curr_anno_dict[(chrom, start)] = anno_info
    return curr_anno_dict


# 前面都统一了注释文件的格式，所以可以直接用
anno_dicts = {
    'effect': get_pos_anno_dict(effect_file),
    'oncokb': get_pos_anno_dict(oncokb_file),
    'cpg1kb': get_pos_anno_dict(cpg1kb_file),
    'cpg2d5kb': get_pos_anno_dict(cpg2d5kb_file),
    'cpg5kb': get_pos_anno_dict(cpg5kb_file),
    'cpg10kb': get_pos_anno_dict(cpg10kb_file),
    'tss1kb': get_pos_anno_dict(tss1kb_file),
    'tss2d5kb': get_pos_anno_dict(tss2d5kb_file),
    'tss5kb': get_pos_anno_dict(tss5kb_file),
    'tss10kb': get_pos_anno_dict(tss10kb_file),
    'rmsk': get_pos_anno_dict(rmsk_file),
    'fullname': get_pos_anno_dict(fullname_file),
    'gc1mb': get_pos_anno_dict(gc1mb_file),
    'tx_exon': get_pos_anno_dict(tx_exon_file),
    'domain': get_pos_anno_dict(domain_file),
}


# 排序一下 Depth
is_cov_dict = {}
with open(is_cov_file) as f:
    next(f)
    for line in f:
        # [240731 MXFA] 新增深度占比和去重 reads 数
        chrom, start, _, umi_num, all_num, all_freq, rmdup_num = line.strip().split('\t')
        is_cov_dict[(chrom, start)] = (int(umi_num), int(all_num), all_freq, int(rmdup_num))
# * 总深度排序 x[1][0]; UMI 排序是 x[1][1]
sorted_is_cov = sorted(is_cov_dict.items(), key=lambda x: x[1][1], reverse=True)


def get_anno(chrom: str, start: str, key: str, idx: int = 0, default: str = '-'):
    """
    获取注释信息. e.g. Intron|FRMD4A
    :param chrom: 染色体
    :param start: 起始位置
    :param key: 注释类型，如 'effect', 'oncokb', 'cpg1kb', 'tss1kb', 'rmsk', 'fullname'
    :param idx: 注释信息的索引，默认为 0
    :param default: 默认值，如果没有找到注释信息，则返回该值
    :return: 注释信息
    """
    val = anno_dicts[key].get((chrom, start))
    if val is None:
        return default
    try:
        return val.split('|')[idx]
    except IndexError:
        return default


# 输出
# [250711] 孟博: 新增 Depth/UMI 比值列
# [250721] 毛博: 新增转录本, 外显子号和基因功能结构域注释
with open(combine_out, 'w') as g:
    headers = ['Chrom', 'Start', 'UMIs', 'Depth', 'Depth/UMI', 'Freq', 'RmdupDepth',
               'Effect', 'Gene', 'FullName', 'Oncogene/TSG',
               'Transcript', 'ExonNum', 'Domain',
               'CpG1KB', 'CpG2.5KB', 'CpG5KB', 'CpG10KB', 'TSS1KB', 'TSS2.5KB', 'TSS5KB', 'TSS10KB',
               'RepName', 'RepClass', 'RepFamily', 'GC1MB']
    g.write('\t'.join(headers) + '\n')
    for (chrom, start), (umi_num, all_num, all_freq, rmdup_num) in sorted_is_cov:
        # ! all_num/umi_num, umi_num 可能为 0 导致除零错误, 因为在前面 gencore 输出的结果 read 信息中 umi 并不在预设 umi 列表中
        all_divide_umi = f'{all_num/umi_num:.4f}' if umi_num != 0 else '0'
        effect = get_anno(chrom, start, 'effect')
        gene = get_anno(chrom, start, 'effect', 1)
        onco = get_anno(chrom, start, 'oncokb')
        fullname = get_anno(chrom, start, 'fullname')
        cpg1kb = anno_dicts['cpg1kb'].get((chrom, start), '-')
        cpg2d5kb = anno_dicts['cpg2d5kb'].get((chrom, start), '-')
        cpg5kb = anno_dicts['cpg5kb'].get((chrom, start), '-')
        cpg10kb = anno_dicts['cpg10kb'].get((chrom, start), '-')
        tss1kb = anno_dicts['tss1kb'].get((chrom, start), '-')
        tss2d5kb = anno_dicts['tss2d5kb'].get((chrom, start), '-')
        tss5kb = anno_dicts['tss5kb'].get((chrom, start), '-')
        tss10kb = anno_dicts['tss10kb'].get((chrom, start), '-')
        gc1mb = anno_dicts['gc1mb'].get((chrom, start), '-')
        rep_name = get_anno(chrom, start, 'rmsk')
        rep_class = get_anno(chrom, start, 'rmsk', 1)
        rep_family = get_anno(chrom, start, 'rmsk', 2)
        transcript = get_anno(chrom, start, 'tx_exon')
        exon_num = get_anno(chrom, start, 'tx_exon', 1)
        domain = get_anno(chrom, start, 'domain')
        g.write('\t'.join([
            chrom, start, str(umi_num), str(all_num), all_divide_umi, all_freq, str(rmdup_num),
            effect, gene, fullname, onco,
            transcript, exon_num, domain,
            cpg1kb, cpg2d5kb, cpg5kb, cpg10kb, tss1kb, tss2d5kb, tss5kb, tss10kb,
            rep_name, rep_class, rep_family, gc1mb
        ]) + '\n')
