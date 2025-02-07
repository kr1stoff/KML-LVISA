"""
过滤 blast 输出符合通过标准的 fastq
   1. Read1 都是 plus (正向)
   2. 3LTR 只允许 1 个 mismatch
   3. 3LTR 引物比对位置过滤. end 位置在 30bp 以内, 引物长度为 25, 允许 5bp 偏移
"""

from dataclasses import dataclass
import sys


@dataclass
class FastqHeader:
    nident: int = 1
    sstrand: int = 3
    slen: int = 7
    qseqid: int = 8
    qend: int = 10


infile = sys.argv[1]
outfile = sys.argv[2]

fh = FastqHeader()
with open(infile, 'r') as f, open(outfile, 'w') as g:
    for line in f:
        lns = line.strip().split('\t')

        # ! 都是正向
        if lns[fh.sstrand] == 'minus':
            continue

        # ! 3LTR mismatch <= 1bp
        if (int(lns[fh.slen]) - int(lns[fh.nident])) > 1:
            continue

        # ! 3LTR 引物比对位置过滤. end 位置在 30bp 以内
        if int(lns[fh.qend]) > 30:
            continue

        g.write(lns[fh.qseqid] + '\n')
