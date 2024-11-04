library(RIdeogram)

# 命令行参数
args <- commandArgs(trailingOnly = TRUE)
combine_file <- args[1] # "anno/SRR17348516.is.combine.tsv"
output_png <- args[2]   # "stats/SRR17348516.chromosome.svg"

load("/data/mengxf/Project/KML240924_lvis_pipeline/script/hg19.cytoBand.overlaid.RData")

# 插入位点统计数据
combine_df <- read.table(combine_file,
                         header = TRUE,
                         sep = "\t")
combine_df <- combine_df[, c(1, 2, 5)]
combine_df$Chrom <- gsub('^chr', '', combine_df$Chrom)
combine_df$End <- as.integer(combine_df$Start + 1)
combine_df$Shape <- 'triangle'
mycolors <- c(
  "E41A1C",
  "377EB8",
  "4DAF4A",
  "984EA3",
  "FF7F00",
  "FFFF33",
  "A65628",
  "F781BF",
  "999999"
)
combine_df$color <- mycolors[match(combine_df$Effect, unique(combine_df$Effect))]
colnames(combine_df) <- c('Chr', 'Start', 'Type', 'End', 'Shape', 'color')
combine_df <- combine_df[, c('Type', 'Shape', 'Chr', 'Start', 'End', 'color')]
# head(combine_df)

# 画图
colorset1 <- c("#FFFFFF",
               "#CCCCCC",
               "#999999",
               "#666666",
               "#333333",
               "#000000")

ideogram(
  karyotype = human_karyotype,
  overlaid = overlaid_tab,
  label = combine_df,
  label_type = "marker",
  colorset1 = colorset1,
  output = output_png
)

convertSVG(
  output_png,
  device = 'png',
  width = 7,
  height = 5,
  dpi = 1000
)
