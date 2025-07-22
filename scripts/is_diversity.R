library(vegan)
library(ineq)


calc_is_diversity <- function(data, output_file) {
    # 计算插入位点多样性的函数
    # 参数:
    #   data: 包含插入位点深度信息的数据框
    #   output_file: 输出文件路径，用于保存多样性统计结果
    # is_abun <- c(10, 5, 2, 8, 6)
    is_abun <- data$Depth
    shannon <- diversity(is_abun, index = "shannon")
    invsimpson <- diversity(is_abun, index = "invsimpson")
    chao1 <- estimateR(is_abun)[2]
    gini_coef <- Gini(is_abun)
    # vegan simpson 是gini-simpson, 经典 simpson = 1 - gini_simpson
    gini_simpson <- diversity(is_abun, index = "simpson")
    stats_df <- data.frame(
        shannon = shannon,
        invsimpson = invsimpson,
        chao1 = chao1,
        gini_coef = gini_coef,
        gini_simpson = gini_simpson
    )
    write.table(
        stats_df,
        file = output_file,
        sep = "\t",
        row.names = F,
        quote = F,
        fileEncoding = "UTF-8"
    )
}


touch_empty_file <- function(output_file) {
    # 创建或清空输出文件的函数
    # 参数:
    #   output_file: 输出文件路径
    if (!file.exists(output_file)) {
        file.create(output_file)
    } else {
        file.remove(output_file)
        file.create(output_file)
    }
}


# I/O
args <- commandArgs(trailingOnly = TRUE)
combine_file <- args[1] # "anno/SRR17348516.is.combine.tsv"
output_file <- args[2] # "stats/SRR17348516.chromosome.svg"

# main
data <- read.table(combine_file, sep = "\t", header = T, quote = "")
if (nrow(data) == 0) {
    touch_empty_file(output_file)
} else {
    calc_is_diversity(data, output_file)
}
