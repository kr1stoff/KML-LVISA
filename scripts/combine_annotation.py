# 合并统计结果到一个表格中

import sys

is_cov_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/is/SRR17348516.is.coverage'
effect_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.effect'
oncokb_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.oncokb'
cpg_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.cpg'
tss_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.tss'
rmsk_file = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.rmsk'
combine_out = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/anno/SRR17348516.is.combine.tsv'


def get_pos_anno_dict(anno_file):
    """
    获取注释信息. #chromosome, start, annotation
    :param anno_file: 注释文件
    :return: 注释字典
    """
    curr_anno_dict = {}
    with open(anno_file) as f:
        for line in f:
            chrom, start, anno_info = line.strip().split('\t')[:3]
            # 重复的位置只要第一个
            if (chrom, start) not in curr_anno_dict:
                curr_anno_dict[(chrom, start)] = anno_info
    return curr_anno_dict


# 前面都统一了注释文件的格式，所以可以直接用
effect_dict = get_pos_anno_dict(effect_file)
oncokb_dict = get_pos_anno_dict(oncokb_file)
cpg_dict = get_pos_anno_dict(cpg_file)
tss_dict = get_pos_anno_dict(tss_file)
rmsk_dict = get_pos_anno_dict(rmsk_file)

# 排序一下 Depth
is_cov_dict = {}
with open(is_cov_file) as f:
    next(f)
    for line in f:
        chrom, start, _, umi_num, all_num = line.strip().split('\t')
        is_cov_dict[(chrom, start)] = (int(umi_num), int(all_num))
# 总深度排序 x[1][0]; UMI 排序是 x[1][1]
sorted_is_cov = sorted(is_cov_dict.items(), key=lambda x: x[1][1], reverse=True)

# 输出
with open(combine_out, 'w') as g:
    headers = ['Chrom', 'Start', 'Depth (UMI)', 'Depth', 'Effect', 'Gene', 'Oncogene/TSG',
               'CpG Name', 'CpG Pos', 'TSS Name', 'TSS Pos', 'Rep Name', 'Rep Class', 'Rep Family']
    g.write('\t'.join(headers) + '\n')

    for (chrom, start), (umi_num, all_num) in sorted_is_cov:
        # effect
        if (chrom, start) in effect_dict:
            effect, gene = effect_dict[(chrom, start)].split('|')
        else:
            effect, gene = '-', '-'
        # oncokb
        if (chrom, start) in oncokb_dict:
            onco, _ = oncokb_dict[(chrom, start)].split('|')
        else:
            onco = '-'
        # cpg
        if (chrom, start) in cpg_dict:
            cpg_name, cpg_pos = cpg_dict[(chrom, start)].split('|')
        else:
            cpg_name, cpg_pos = '-', '-'
        # tss
        if (chrom, start) in tss_dict:
            tss_name, tss_pos = tss_dict[(chrom, start)].split('|')
        else:
            tss_name, tss_pos = '-', '-'
        # rmsk
        # repName, repClass, repFamily
        if (chrom, start) in rmsk_dict:
            rep_name, rep_class, rep_family = rmsk_dict[(chrom, start)].split('|')
        else:
            rep_name, rep_class, rep_family = '-', '-', '-'

        g.write('\t'.join([chrom, start, str(umi_num), str(all_num), effect, gene, onco, cpg_name,
                cpg_pos, tss_name, tss_pos, rep_name, rep_class, rep_family]) + '\n')
