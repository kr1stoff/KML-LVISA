library("CMplot")

# I/O
args <- commandArgs(trailingOnly = TRUE)
combine_file <- args[1] # "anno/SRR17348516.is.combine.tsv"
output_png <- args[2] # "stats/SRR17348516.chromosome_density.png"

data_comb <- read.table(combine_file, header = TRUE, sep = "\t")[, 1:4]
data_comb <- cbind(paste("SNP", rownames(data_comb), sep = "_"), data_comb)
colnames(data_comb) <- c("SNP", "Chromosome", "Position", "trait1", "trait2")

png(output_png, units = "in", res = 300, width = 9, height = 6)
CMplot(data_comb,
    plot.type = "d", bin.size = 1e6, chr.den.col = c("darkgreen", "yellow", "red"),
    file = "png", file.name = "", dpi = 300, main = "Chromosom Integration Site Density",
    file.output = FALSE, verbose = TRUE, width = 9, height = 6
)
dev.off()
