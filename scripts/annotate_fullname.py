# 根式 HGNC 注释基因的全名

import sys

snpeff_parse = sys.argv[1]
fullname_out = sys.argv[2]
hgnc_db = sys.argv[3]  # /data/mengxf/Database/genome/hg19/HGNC_gene_with_protein_product.tsv

symbol_name_dict = {}

with open(hgnc_db) as f:
    header = next(f)
    hd = {h[1]: h[0] for h in enumerate(header.strip().split('\t'))}
    for line in f:
        lns = line.split('\t')
        symbol = lns[hd['symbol']]
        name = lns[hd['name']]
        symbol_name_dict[symbol] = name
        # * alias_symbol
        alias_symbols = lns[hd['alias_symbol']].split('|')
        for alias_symbol in alias_symbols:
            symbol_name_dict[alias_symbol] = name
        # * prev_symbol
        prev_symbols = lns[hd['prev_symbol']].split('|')
        for prev_symbol in prev_symbols:
            symbol_name_dict[prev_symbol] = name

with open(snpeff_parse) as f, open(fullname_out, 'w') as g:
    for line in f:
        chrom, start, effect_and_gene = line.strip().split('\t')
        _, gene = effect_and_gene.split('|')

        if gene in symbol_name_dict:
            g.write(f'{chrom}\t{start}\t{symbol_name_dict[gene]}|{gene}\n')
