import re
from random import sample
import sys

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input[0]) as f, open(snakemake.output[0], 'w') as g:
    for line in f:
        chrom, start, end = re.findall(r'>(chr.*?):(\d+)-(\d+)', line)[0]
        refbase = next(f).strip().upper()
        altbase = sample('ATGC'.replace(refbase, ''), 1)[0]
        # ! BED: 0-based; VCF: 1-based
        print(chrom, str(int(start)+1), end, refbase, altbase, '-', '-', '-', file=g, sep='\t')
