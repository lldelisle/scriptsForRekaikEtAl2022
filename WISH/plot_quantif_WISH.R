# Load packages
if (!"reshape" %in% installed.packages()) {
  install.packages("reshape", repos = "https://stat.ethz.ch/CRAN/")
}
library(reshape)
# if (!"plyr" %in% installed.packages()) {
#   install.packages("plyr", repos = "https://stat.ethz.ch/CRAN/")
# }
# library(plyr)
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
quantifs.file <- "WISH/quantification.txt"
dir.plot <- dirname(quantifs.file)

quantifs.table <- read.delim(quantifs.file, skip = 3, header = F)
genotypes <- c("wt", "Del(CBS1-2)")
colors <- c("#3e8c41", "#ed1d24")
names(colors) <- genotypes
genes <- c("Hoxd4", "Hoxd9")
colnames(quantifs.table) <- paste0(rep(c("domain.size", "stembryo.length", "ratio"), length(genotypes) * length(genes)),
                                   "_",
                                   rep(rep(genotypes, each = 3), length(genes)),
                                   "_",
                                   rep(genes, each = 3 * length(genotypes)))
for(gene in genes) {
  for(genotype in genotypes) {
    if (! all(round(quantifs.table[, paste0("ratio_", genotype, "_", gene)], 2)
              == round(quantifs.table[, paste0("domain.size_", genotype, "_", gene)] /
                       quantifs.table[, paste0("stembryo.length_", genotype, "_", gene)], 2), na.rm = T) ){
      stop(gene, genotype, "ratio is incorrect")
    }
  }
}
quantifs.data <- melt(quantifs.table[, grep("ratio", colnames(quantifs.table))])
quantifs.data$genotype <- sapply(strsplit(as.character(quantifs.data$variable), "_"), "[[", 2)
quantifs.data$genotype <- factor(quantifs.data$genotype, levels = genotypes)
quantifs.data$gene <- sapply(strsplit(as.character(quantifs.data$variable), "_"), "[[", 3)
quantifs.data$gene <- factor(quantifs.data$gene, levels = rev(genes))

simple.fun <- function(v){return(list(ymin = mean(v) - sd(v),
                                      ymax = mean(v) + sd(v)))}

ggplot(quantifs.data, aes(x = genotype, y = value, color = genotype)) +
  geom_jitter(width = 0.3) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.5)  +
  stat_summary(fun.min = mean, fun.max = mean,
               geom="errorbar", color="black", width=0.8)  +
  facet_grid(gene ~ .) +
  theme_classic() +
  ylim(0, 1) +
  scale_color_manual(values = colors)
ggsave(file.path(dir.plot, "FigS13.pdf"), height = 7, width = 4)
