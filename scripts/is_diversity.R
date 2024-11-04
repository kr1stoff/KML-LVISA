library(vegan)
library(ineq)

# I/O
args <- commandArgs(trailingOnly = TRUE)
combine_file <- args[1] # "anno/SRR17348516.is.combine.tsv"
output_file <- args[2]   # "stats/SRR17348516.chromosome.svg"

# is_abun <- c(10, 5, 2, 8, 6)
data <- read.table(combine_file, sep = '\t', header = T)
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
  sep = '\t',
  row.names = F,
  quote = F,
  fileEncoding = "UTF-8"
)
