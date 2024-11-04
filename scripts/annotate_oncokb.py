# 根据 OncoKB 注释基因列表

import sys

oncokb_db = '/data/mengxf/Database/genome/hg19/OncoKBcancerGeneList.tsv'
snpeff_parse = sys.argv[1]
oncokb_out = sys.argv[2]

# 获取 oncogen 和 tsg 列表
onco_list, tsg_lst = [], []
with open(oncokb_db) as f:
    header = next(f)
    hd = {h[1]: h[0] for h in enumerate(header.strip().split('\t'))}

    for line in f:
        lns = line.split('\t')
        gene = lns[hd['Hugo Symbol']]

        if lns[hd['Is Oncogene']] == 'Yes':
            onco_list.append(gene)
        if lns[hd['Is Tumor Suppressor Gene']] == 'Yes':
            tsg_lst.append(gene)

        # snpEff 没有用到基因别名
        # alias = lns[hd['Gene Aliases']].strip().split(', ')
        # genes = [gene] + alias
        # for gn in genes:
        #     if lns[hd['Is Oncogene']] == 'Yes':
        #         onco_list.append(gn)
        #     if lns[hd['Is Tumor Suppressor Gene']] == 'Yes':
        #         tsg_lst.append(gn)


with open(snpeff_parse) as f, open(oncokb_out, 'w') as g:
    for line in f:
        chrom, start, effect_and_gene = line.strip().split('\t')
        gene = effect_and_gene.split('|')[1]
        if gene in onco_list:
            g.write(f'{chrom}\t{start}\toncogene|{gene}\n')
        elif gene in tsg_lst:
            g.write(f'{chrom}\t{start}\ttsg|{gene}\n')
