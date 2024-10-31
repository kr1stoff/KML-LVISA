# 过滤 blast 输出符合通过标准的 fastq header
#   1. Read1 都是 plus (正向)
#   2. 3LTR 只允许 1 个 mismatch

from dataclasses import dataclass
import sys


@dataclass
class FastqHeader:
    nident: int = 1
    sstrand: int = 3
    slen: int = 7
    qseqid: int = 8


infile = sys.argv[1]  # '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/3ltr/read1.blastout '
outfile = sys.argv[2]  # '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/3ltr/qseqid_read1.txt'

fh = FastqHeader()
with open(infile, 'r') as f, open(outfile, 'w') as g:
    for line in f:
        lns = line.strip().split('\t')

        # 都是正向
        if lns[fh.sstrand] == 'minus':
            continue

        # 3LTR mismatch <= 1bp
        if (int(lns[fh.slen]) - int(lns[fh.nident])) > 1:
            continue

        g.write(lns[fh.qseqid] + '\n')
