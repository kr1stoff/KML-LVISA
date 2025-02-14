library("CMplot")


draw_cmplot <- function(data_comb, output_png) {
    # 计算插入位点密度的函数
    # 参数:
    #   data: 包含插入位点信息的数据框
    #   output_png: 输出文件路径，用于保存密度图
    # 计算插入位点密度
    data_comb <- cbind(paste("SNP", rownames(data_comb), sep = "_"), data_comb)
    colnames(data_comb) <- c("SNP", "Chromosome", "Position", "trait1", "trait2")

    png(output_png, units = "in", res = 300, width = 9, height = 6)
    CMplot(data_comb,
        plot.type = "d", bin.size = 1e6, chr.den.col = c("darkgreen", "yellow", "red"),
        file = "png", file.name = "", dpi = 300, main = "Chromosom Integration Site Density",
        file.output = FALSE, verbose = TRUE, width = 9, height = 6
    )
    dev.off()
}

# 创建或清空输出文件的函数
# 参数:
#   output_file: 输出文件路径
touch_empty_file <- function(output_file) {
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
output_png <- args[2] # "stats/SRR17348516.chromosome_density.png"

# main
data_comb <- read.table(combine_file, header = TRUE, sep = "\t")[, 1:4]
# ! NTC 位点数太少会报错，改成 try-catch
tryCatch(
    {
        draw_cmplot(data_comb, output_png)
    },
    error = function(e) {
        message("An error occurred: ", e$message)
        touch_empty_file(output_png)
    }
)
