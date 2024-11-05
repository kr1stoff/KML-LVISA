# 数据质控
fastqc -t 4 --extract -o fastqc ../../rawdata/SRR17348516_1.fastq.gz ../../rawdata/SRR17348516_2.fastq.gz
# UMI, R2测穿反向3LTR
fastp -w 8 -q 15 -u 40 -l 82 \
    --adapter_sequence_r2 AAATGGTCTGAGGGATCTCTAGT \
    --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction \
    --umi --umi_loc=read2 --umi_len=10 \
    -i ../../rawdata/SRR17348516_1.fastq.gz \
    -I ../../rawdata/SRR17348516_2.fastq.gz \
    -o fastp/SRR17348516.clean.R1.fq \
    -O fastp/SRR17348516.clean.R2.fq \
    -j fastp/SRR17348516.json \
    -h fastp/SRR17348516.html

# 3LTR
seqtk seq -A fastp/SRR17348516.clean.R1.fq >3ltr/SRR17348516.clean.R1.fa
blastn -task blastn-short \
    -word_size 7 \
    -query 3ltr/SRR17348516.clean.R1.fa \
    -db /data/mengxf/GitHub/KML-LVISA/assets/3LTR/3LTR.fa \
    -out 3ltr/blast.out \
    -num_threads 16 \
    -max_hsps 100 \
    -max_target_seqs 1000 \
    -outfmt "6 length nident pident sstrand sseqid sstart send slen qseqid qstart qend qlen qcovs evalue bitscore"
python /data/mengxf/GitHub/KML-LVISA/wf-lvisa/scripts/get_fastq_header_from_blast.py
seqkit grep -f 3ltr/qseqid.txt fastp/SRR17348516.clean.R1.fq -o 3ltr/SRR17348516.1.3ltr.fq
seqkit grep -f 3ltr/qseqid.txt fastp/SRR17348516.clean.R2.fq -o 3ltr/SRR17348516.2.3ltr.fq

# 比对
bwa mem -t 16 -M -Y -R '@RG\tID:SRR17348516\tSM:SRR17348516' \
    /data/mengxf/Database/genome/hg19/hg19.fa \
    3ltr/SRR17348516.1.3ltr.fq \
    3ltr/SRR17348516.2.3ltr.fq |
    samtools view -@ 4 -hbS - |
    samtools sort -@ 4 -o align/SRR17348516.bam -
samtools stat align/SRR17348516.bam | grep ^SN | cut -f 2- >align/SRR17348516.bam.stat
# 过滤 1.mapq < 30;  2.map length < 55;  3. flag 4 8 256 512 2048
samtools view -@ 4 -hbS -q 30 -m 55 -F 2828 align/SRR17348516.bam -o align/SRR17348516.filter.bam
samtools stat align/SRR17348516.filter.bam | grep ^SN | cut -f 2- >align/SRR17348516.filter.bam.stat

# umi  gencore 貌似只能去重不能统计
gencore -i align/SRR17348516.filter.bam \
    -o umi/SRR17348516.umi.bam \
    -r /data/mengxf/Database/genome/hg19/hg19.fa \
    -h umi/SRR17348516.umi.html \
    -j umi/SRR17348516.umi.json
samtools stat umi/SRR17348516.umi.bam | grep ^SN | cut -f 2- >umi/SRR17348516.umi.bam.stat
# 过滤小于3条 reads 的 UMI
/home/mengxf/miniforge3/envs/python3.8/bin/python \
    /data/mengxf/Project/KML240924_lvis_pipeline/script/filter_umi_bam.py \
    umi/SRR17348516.umi.bam \
    umi/SRR17348516.umi.filter.bam

# IS 计数
# umi bed 区域
bedtools bamtobed -i umi/SRR17348516.umi.filter.bam | bedtools merge >is/SRR17348516.umi.bed
# uniq umi bam 深度. coverage 解释: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html#default-behavior
bedtools coverage -a is/SRR17348516.umi.bed -b umi/SRR17348516.umi.filter.bam >is/SRR17348516.umi.coverage
# 总 bam 深度
bedtools coverage -a is/SRR17348516.umi.bed -b align/SRR17348516.filter.bam >is/SRR17348516.all.coverage
# 过滤 1个loci < 2 UMI，support reads < 10
/home/mengxf/miniforge3/envs/python3.8/bin/python \
    /data/mengxf/Project/KML240924_lvis_pipeline/script/filter_and_stat_is_by_coverage.py \
    is/SRR17348516.umi.coverage \
    is/SRR17348516.all.coverage \
    is/SRR17348516.is.coverage

# 注释
# IS 只要插入位置那 1bp 的碱基
awk -F '\t' -v OFS='\t' '{print $1,$2,$2+1}' is/SRR17348516.is.coverage >anno/SRR17348516.is.bed
snpEff -dataDir /data/mengxf/Database/snpEff \
    -i bed -chr chr -geneId -canon -noStats \
    hg19 anno/SRR17348516.is.bed >anno/SRR17348516.is.bed.snpEff
# 注释 position, gene, effect
/home/mengxf/miniforge3/envs/python3.8/bin/python \
    /data/mengxf/Project/KML240924_lvis_pipeline/script/annotate_effect.py \
    anno/SRR17348516.is.bed.snpEff \
    anno/SRR17348516.is.effect
# 注释 oncokb
/home/mengxf/miniforge3/envs/python3.8/bin/python \
    /data/mengxf/Project/KML240924_lvis_pipeline/script/annotate_oncokb.py \
    anno/SRR17348516.is.bed.snpEff \
    anno/SRR17348516.is.effect
# CpG 要去重, CpG(10kb) 之间有 overlap. 后面都做去重
bedtools intersect -wb \
    -a anno/SRR17348516.is.bed \
    -b /data/mengxf/Database/genome/hg19/hg19.cpg10kb.bed |
    cut -f1,2,7 >anno/SRR17348516.is.cpg
# TSS
bedtools intersect -wb \
    -a anno/SRR17348516.is.bed \
    -b /data/mengxf/Database/genome/hg19/hg19.switchDbTss10kb.bed |
    cut -f1,2,7 >anno/SRR17348516.is.tss
# Repeat
# repName, repClass, repFamily
bedtools intersect -wb \
    -a anno/SRR17348516.is.bed \
    -b /data/mengxf/Database/genome/hg19/hg19.rmsk.bed |
    cut -f1,2,7 >anno/SRR17348516.is.rmsk
# 合并注释结果
/home/mengxf/miniforge3/envs/python3.8/bin/python /data/mengxf/Project/KML240924_lvis_pipeline/script/combine_annotation.py

# 统计
# 染色体分布图

# 多样性
/home/mengxf/miniforge3/envs/R4.3.1/bin/Rscript \
    /data/mengxf/Project/KML240924_lvis_pipeline/script/is_diversity.R \
    anno/SRR17348516.is.combine.tsv \
    stats/is.diversity.txt
# Top 10 IS 位点柱状图
/home/mengxf/miniforge3/envs/python3.8/bin/python /data/mengxf/Project/KML240924_lvis_pipeline/script/top10_is.py
# 染色体插入位点占比柱状图 (均一化化染色体长度)
/home/mengxf/miniforge3/envs/python3.8/bin/python /data/mengxf/Project/KML240924_lvis_pipeline/script/chromosome_distribute_barplot.py
# snpEff effect 饼图
/home/mengxf/miniforge3/envs/python3.8/bin/python /data/mengxf/Project/KML240924_lvis_pipeline/script/effect_pie.py
# Repeat repClass 饼图
/home/mengxf/miniforge3/envs/python3.8/bin/python /data/mengxf/Project/KML240924_lvis_pipeline/script/reclass_pie.py
