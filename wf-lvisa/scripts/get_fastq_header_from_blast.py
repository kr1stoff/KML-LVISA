from dataclasses import dataclass


@dataclass
class FastqHeader:
    nident: int = 1
    sstrand: int = 3
    slen: int = 7
    qseqid: int = 8


infile = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/3ltr/blast.out'
outfile = '/data/mengxf/Project/KML240924_lvis_pipeline/result/240929/3ltr/qseqid.txt'

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
