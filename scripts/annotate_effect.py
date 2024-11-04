# 使用 snpEff 信息注释插入位点
# chromosome, start, effect, gene
import sys

snpeff_out = sys.argv[1]
anno_out = sys.argv[2]


def parse_anno_info(anno: str):
    # Name;Effect|Gene|BioType
    egb = anno.split(';')[1]
    eff = egb.split('|')[0].split(':')[0]
    gene = '-'

    # 根据不同 effect 类型，提取 gene 和 biotype 信息
    # Intergenic:ZNF536-TSHZ3
    # Gene:TAF1:pseudogene
    if eff in ['Intergenic', 'Gene']:
        gene = egb.split(':')[1]
    # Intron:2:5:RETAINED-RETAINED|Transcript:NM_016205.3|Gene:PDGFC:pseudogene
    # Downstream:4810|Transcript:NM_198542.3|Gene:ZNF773:protein_coding;Gene:ZNF773:protein_coding
    # Exon:9:9:ALTTENATIVE_POLY_A|Transcript:NM_001013437.2|Gene:SEH1L:protein_coding
    # Upstream:671|Transcript:NM_145865.3|Gene:ANKS4B:protein_coding
    elif eff in ['Intron', 'Downstream', 'Exon', 'Upstream']:
        _, _, gene_raw = egb.split('|')
        gene = gene_raw.split(':')[1]

    return eff, gene


with open(snpeff_out) as f, open(anno_out, 'w') as g:
    for line in f:
        if line.startswith('#'):
            continue

        lns = line.strip().split('\t')
        chrom, start, _, anno = lns
        eff, gene = parse_anno_info(anno)
        g.write(f'{chrom}\t{start}\t{eff}|{gene}\n')
