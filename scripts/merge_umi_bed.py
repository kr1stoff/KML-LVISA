import sys
import json
import Levenshtein


def assign_umi_seq_to_id(umi_seq):
    """根据 UMI 序列计算 UMI ID"""
    for k, v in umi_dict.items():
        if Levenshtein.distance(umi_seq, k) <= 1:
            return v
    return 'UMI0'


# options
strand_bed = sys.argv[1]
qname_bed = sys.argv[2]
umi_json = sys.argv[3]
out_bed = sys.argv[4]

# main
umi_dict = json.load(open(umi_json, 'r'))

with open(strand_bed) as f, open(qname_bed) as f2, open(out_bed, 'w') as g:
    for sline in f:
        qline = next(f2)
        slns = sline.strip().split('\t')
        qlns = qline.strip().split('\t')
        # 位置
        spos = slns[:3]
        qpos = qlns[:3]
        if qpos != qpos:
            raise Exception('strand 和 qname bed 位置不一致!')
        # 链
        strand = slns[3]
        # UMI
        qname = qlns[3]
        raw_umi_seqs = [q.split('/')[0].split(':')[-1] for q in qname.strip().split(',')]
        umis = list(set([assign_umi_seq_to_id(umi_seq) for umi_seq in raw_umi_seqs]))
        umi_num = len(list(filter(lambda x: x != 'UMI0', umis)))
        # 输出
        g.write('\t'.join(spos + [strand, str(umi_num), ','.join(umis)]) + '\n')
