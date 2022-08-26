# Load packages
if (!"reshape" %in% installed.packages()) {
  install.packages("reshape", repos = "https://stat.ethz.ch/CRAN/")
}
library(reshape)
if (!"plyr" %in% installed.packages()) {
  install.packages("plyr", repos = "https://stat.ethz.ch/CRAN/")
}
library(plyr)
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
if (!"ggh4x" %in% installed.packages()) {
  install.packages("ggh4x", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggh4x)

# Set paths
directory.with.FPKM <- commandArgs(TRUE)[1]
directory.with.plots <- commandArgs(TRUE)[2]
directory.with.data <- file.path(directory.with.plots, "..", "raw_values")
# Create non-existing output directories
if (!dir.exists(directory.with.plots)) {
  dir.create(directory.with.plots, recursive = T, showWarnings = F)
}
if (!dir.exists(directory.with.data)) {
  dir.create(directory.with.data, recursive = T, showWarnings = F)
}

# Function to plot error bars
simple.fun <- function(v){return(list(ymin = mean(v) - sd(v),
                                      ymax = mean(v) + sd(v)))}

# Generate table with all FPKMs:
all.FPKM.files <- list.files(path = directory.with.FPKM, pattern = "FPKM.txt.gz", recursive = T, full.names = T)
all.FPKM.values <- do.call(rbind, lapply(all.FPKM.files, function(fn){
  temp.df <- read.delim(fn)
  temp.df <- subset(temp.df, select = c("gene_id", "gene_short_name", "FPKM"))
  temp.df$sample <- gsub("_FPKM.txt.gz", "", basename(fn))
  return(temp.df)
}))
# Write all FPKMs in table (summing for identical gene_ids)
write.table(cast(all.FPKM.values, gene_id + gene_short_name ~ sample, value = "FPKM", fun.aggregate = sum),
            file.path(directory.with.data, paste0("all_FPKM.txt")),
            sep = "\t", quote = F, row.names = F)
system(paste("gzip -f", file.path(directory.with.data, paste0("all_FPKM.txt"))))


# Figure 1E:
# Get FPKM values for heatmap:
all.genes.needed <- paste0("Hoxd", 1:13)
all.FPKM.files <- grep("reptc", list.files(path = directory.with.FPKM, pattern = "FPKM.txt.gz", recursive = T, full.names = T),
                       value = T)
all.FPKM.values <- do.call(rbind, lapply(all.FPKM.files, function(fn){
  temp.df <- read.delim(fn)
  temp.df <- subset(temp.df, subset = gene_short_name %in% all.genes.needed,
                    select = c("gene_short_name", "FPKM"))
  temp.df$sample <- gsub("_FPKM.txt.gz", "", basename(fn))
  return(temp.df)
}))

# Extract the samplesplan from names
samplesplan <- data.frame(sample = unique(all.FPKM.values$sample))
samplesplan$time <- sapply(strsplit(samplesplan$sample, "_"), "[[", 2)

# Plot the heatmap (Fig1E):
summary.df <- ddply(subset(merge(all.FPKM.values, samplesplan), grepl("Hoxd", gene_short_name)), .(gene_short_name, time), summarize,
                    mean_value = mean(FPKM))
summary.df$time <- factor(summary.df$time, levels = paste0(seq(168, 72, -12), "h"))
summary.df$gene_short_name <- factor(gsub("Hox", "", summary.df$gene_short_name),
                                     levels = paste0("d", 13:1))
ggplot(summary.df, aes(gene_short_name, time, fill = mean_value)) +
  geom_tile() +
  scale_fill_distiller("FPKM",
                       palette = "BuPu",
                       direction = 0) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "italic")) +
  ggtitle("RNA-seq")

ggsave(file.path(directory.with.plots, "Fig1E.pdf"), height = 3, width = 4)
# Export data:
write.table(cast(summary.df, time ~ gene_short_name, value = "mean_value"),
            file.path(directory.with.data, paste0("Fig1E_raw_data.txt")),
            sep = "\t", quote = F, row.names = F)

# Figure S2Cb:
ratios <- data.frame(what = "FPKM", time = c("96h", "120h"),
                     ratio = c(summary.df$mean_value[summary.df$time == "96h" & summary.df$gene_short_name == "d9"] / 
                                 summary.df$mean_value[summary.df$time == "96h" & summary.df$gene_short_name == "d4"],
                               summary.df$mean_value[summary.df$time == "120h" & summary.df$gene_short_name == "d9"] / 
                                 summary.df$mean_value[summary.df$time == "120h" & summary.df$gene_short_name == "d4"]))
ratios$time <- factor(ratios$time, levels = c("96h", "120h"))
ggplot(ratios, aes(x = time, y = ratio)) +
  geom_line(aes(group = what)) +
  facet_grid(. ~ what) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("") +
  ylab("d9/d4 FPKM ratio") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))
ggsave(file.path(directory.with.plots, paste0("FigS2Cb.pdf")), height = 3, width = 2)

# Figure 7C
# Get FPKM values for heatmap Del(CBS1-5):
all.genes.needed <- paste0("Hoxd", 1:13)
all.FPKM.files <- list.files(path = directory.with.FPKM, pattern = "Del\\(CBS1-5\\).*FPKM.txt.gz", recursive = T, full.names = T)
all.FPKM.values <- do.call(rbind, lapply(all.FPKM.files, function(fn){
  temp.df <- read.delim(fn)
  temp.df <- subset(temp.df, subset = gene_short_name %in% all.genes.needed,
                    select = c("gene_short_name", "FPKM"))
  temp.df$sample <- gsub("_FPKM.txt.gz", "", basename(fn))
  return(temp.df)
}))

samplesplan <- data.frame(sample = unique(all.FPKM.values$sample))
samplesplan$time <- sapply(strsplit(samplesplan$sample, "_"), "[[", 2)

# Plot the heatmap:
summary.df <- ddply(merge(all.FPKM.values, samplesplan), .(gene_short_name, time), summarize,
                    mean_value = mean(FPKM))
summary.df$time <- factor(summary.df$time, levels = paste0(seq(168, 72, -12), "h"))
summary.df$gene_short_name <- factor(gsub("Hox", "", summary.df$gene_short_name),
                                     levels = paste0("d", 13:1))
ggplot(summary.df, aes(gene_short_name, time, fill = mean_value)) +
  geom_tile() +
  scale_fill_distiller("FPKM",
                       palette = "BuPu",
                       direction = 0) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "italic")) +
  ggtitle("Del(CBS1-5)")

ggsave(file.path(directory.with.plots, "Fig7C.pdf"), height = 2, width = 4)
# Export data:
write.table(cast(summary.df, time ~ gene_short_name, value = "mean_value"),
            file.path(directory.with.data, paste0("Fig7C_raw_data.txt")),
            sep = "\t", quote = F, row.names = F)


# Normalization matrix
# The FPKM of the gene used as name is
# normalized by the mean of the genes
# in this list
normalization.matrix <- list("d1" = paste0("a", 1:2),
                             "d3" = "a3",
                             "d4" = c("a5", "c4"),
                             "d8" = "c8",
                             "d9" = paste0(c("a", "c"), 9),
                             "d10" = paste0(c("a", "c"), 10),
                             "d11" = paste0(c("a", "c"), 11),
                             "a1" = "b1",
                             "a5" = paste0(c("b", "c"), 5),
                             "a9" = "c9")

gene.colors <- c('d13' = "#752265", 'd11' = "#b02026", 'd10' = "#087b56",
                 'd9' = "#3c58a7", "d8" = "#ee3124", "d4" = "#32b44a", "d3" = "#9753a1",
                 "d1" = "#f7931d", "a9" = "#040707", "a5" = "#a5782b", "a1" = "#293386")
gene.shapes <- c('d13' = 8, 'd11' = 6, 'd10' = 5, "d9" = 16, "d8" = 15, "d4" = 17, "d3" = 17, "d1" = 18, "a9" = 19, "a5" = 15, "a1" = 17)


# Store all necessary genes in an array
all.genes.needed <- c(paste0("Hox", names(normalization.matrix)),
                      paste0("Hox", unlist(normalization.matrix)),
                      "Hoxd13")
# Get all FPKM files which requires normalization
all.FPKM.files <- grep("reptc", list.files(path = directory.with.FPKM, pattern = "FPKM.txt.gz", recursive = T, full.names = T),
                       invert = T, value = T)
# Get the FPKM values for these genes
all.FPKM.values <- do.call(rbind, lapply(all.FPKM.files, function(fn){
  temp.df <- read.delim(fn)
  temp.df <- subset(temp.df, subset = gene_short_name %in% all.genes.needed,
                    select = c("gene_short_name", "FPKM"))
  temp.df$sample <- gsub("_FPKM.txt.gz", "", basename(fn))
  return(temp.df)
}))
# Process sample names to get info
samplesplan <- data.frame(sample = unique(all.FPKM.values$sample))
samplesplan$genotype <- sapply(strsplit(samplesplan$sample, "_"), "[[", 1)
samplesplan$time <- sapply(strsplit(samplesplan$sample, "_"), "[[", 2)
samplesplan$rep <- sapply(strsplit(samplesplan$sample, "_"), "[[", 4)

# Switch to samples as row and genes as columns:
all.FPKM.values.v1 <- cast(all.FPKM.values, sample ~ gene_short_name, value = "FPKM")
# Compute the normalization and store it in 'gene_ratio'
for (gene in names(normalization.matrix)) {
  if (length(normalization.matrix[[gene]]) == 1) {
    denom <- all.FPKM.values.v1[, paste0("Hox", normalization.matrix[[gene]])]
  } else {
    denom <- apply(all.FPKM.values.v1[, paste0("Hox", normalization.matrix[[gene]])], 1, mean)
  }
  all.FPKM.values.v1[, paste0("Hox", gene, "_ratio")] <- all.FPKM.values.v1[, paste0("Hox", gene)] / denom
}
# Switch back to one value per line:
all.FPKM.values.v2 <- melt(as.data.frame(all.FPKM.values.v1), id.vars = "sample", measure.vars = setdiff(colnames(all.FPKM.values.v1), "sample"), variable_name = "gene")
# Add the sample info
all.FPKM.values.v2 <- merge(all.FPKM.values.v2, samplesplan)

# Normalize by the corresponding wt replicates:
# For Del(sub-TAD1) and Del(CBS1-5), the best controls are wt rep3 and rep4
# However, they are not present at 144h
# Unfortunately Hoxa3 is to 0 in wt_120h_RNAseq_rep4
# Thus only rep3 is used in this case.
all.FPKM.values.v2 <- rbind(
  ddply(subset(all.FPKM.values.v2, !rep %in% c("rep3", "rep4") & (!genotype %in% c("Del(sub-TAD1)", "Del(CBS1-5)") | time == "144h")), .(gene, time), mutate,
        control2rep = mean(value[rep %in% c("rep1", "rep2") & genotype == "wt"])
  ),
  ddply(subset(all.FPKM.values.v2, (gene != "Hoxd3_ratio" | time != "120h") & (rep %in% c("rep3", "rep4") | (genotype %in% c("Del(sub-TAD1)", "Del(CBS1-5)") & time != "144h"))), .(gene, time), mutate,
        control2rep = mean(value[rep %in% c("rep3", "rep4") & genotype == "wt"])
  ),
  # Correct the control2rep for Hoxd3_ratio:
  ddply(subset(all.FPKM.values.v2, gene == "Hoxd3_ratio" & time == "120h" & (rep %in% c("rep3", "rep4") | (genotype %in% c("Del(sub-TAD1)", "Del(CBS1-5)") & time != "144h"))), .(gene, time), mutate,
        control2rep = mean(value[rep %in% c("rep3") & genotype == "wt"])
  )
)

# Compute the ratio with the wt samples
all.FPKM.values.v2 <- ddply(all.FPKM.values.v2, .(sample, gene), mutate,
                            ratio_wt = value / control2rep)

to.plot.1 <- subset(all.FPKM.values.v2, gene %in% paste0("Hox", c("d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"), "_ratio") &
                      !genotype %in% c("Del(CBS4)", "Del(CBS1-5)") & time != "144h")
to.plot.1$gene <- gsub("Hox", "", gsub("_ratio", "", to.plot.1$gene))
to.plot.1$gene <- factor(to.plot.1$gene, levels = c("d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.1$time <- factor(to.plot.1$time, levels = c("96h", "120h", "144h"))

figures <- list("Fig5B" = list(geno = "Del(CBS1)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2, 4),
                                                                limits = c(0, 4),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)))),
                "Fig5D" = list(geno = "Del(CBS1-2)",
                               scales = list(scale_y_continuous(breaks = c(2, 4, 6, 8),
                                                                limits = c(2, 8),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)))),
                "Fig5F" = list(geno = "Del(CBS2)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2, 4),
                                                                limits = c(0, 4),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)))),
                "Fig5H" = list(geno = "Ins(2xCBS-d4d8)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)))),
                "Fig6B" = list(geno = "Del(d1-d4)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)))),
                "Fig6E" = list(geno = "Del(sub-TAD1)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0))))
)

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.1, genotype %in% c("wt", figures[[fig.name]][["geno"]]) & is.finite(ratio_wt))
  if (figures[[fig.name]][["geno"]] != "Del(sub-TAD1)") {
    df.stats <- subset(df.stats, rep %in% c("rep1", "rep2"))
  } else {
    df.stats <- subset(df.stats, genotype == "Del(sub-TAD1)" |
                         rep %in% c("rep3", "rep4"))
  }
  if (fig.name == "Fig5D") {
    df.stats$time <- as.character(df.stats$time)
    temp.df <- subset(df.stats, time == "96h")
    df.stats$time[df.stats$time == "96h"] <- "high 96h"
    temp.df$time <- "low 96h"
    df.stats <- rbind(df.stats, temp.df)
    df.stats$time <- factor(df.stats$time, levels = c("high 96h", "low 96h", "120h", "144h"))
  }
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  ns <- table(df.stats$split)
  df.stats.split <- split(subset(df.stats, split %in% names(ns[ns > 3])), factor(subset(df.stats, split %in% names(ns[ns > 3]))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("Normalized FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, ".pdf")), g, height = 5, width = 7)
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}


to.plot.1.noratio <- subset(all.FPKM.values.v2, gene %in% paste0("Hox", c("d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1")) &
                              !genotype %in% c("Del(CBS4)", "Del(CBS1-5)") &  time != "144h")
to.plot.1.noratio$gene <- gsub("Hox", "", to.plot.1.noratio$gene)
to.plot.1.noratio$gene <- factor(to.plot.1.noratio$gene, levels = c("d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.1.noratio$time <- factor(to.plot.1.noratio$time, levels = c("96h", "120h", "144h"))

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.1.noratio, genotype %in% c("wt", figures[[fig.name]][["geno"]]) & is.finite(ratio_wt))
  if (figures[[fig.name]][["geno"]] != "Del(sub-TAD1)") {
    df.stats <- subset(df.stats, rep %in% c("rep1", "rep2"))
  } else {
    df.stats <- subset(df.stats, genotype == "Del(sub-TAD1)" |
                         rep %in% c("rep3", "rep4"))
  }
  if (fig.name == "Fig5D") {
    df.stats$time <- as.character(df.stats$time)
    temp.df <- subset(df.stats, time == "96h")
    df.stats$time[df.stats$time == "96h"] <- "high 96h"
    temp.df$time <- "low 96h"
    df.stats <- rbind(df.stats, temp.df)
    df.stats$time <- factor(df.stats$time, levels = c("high 96h", "low 96h", "120h", "144h"))
  }
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  ns <- table(df.stats$split)
  df.stats.split <- split(subset(df.stats, split %in% names(ns[ns > 3])), factor(subset(df.stats, split %in% names(ns[ns > 3]))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, "_not_norm.pdf")), g, height = 5, width = 7)
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "not_norm_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}

to.plot.2 <- subset(all.FPKM.values.v2, 
                    gene %in% c("Hoxd13", paste0("Hox", c("d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"), "_ratio")) &
                      (genotype == "Del(CBS1-5)"  | (genotype == "wt" & rep %in% c("rep1", "rep2") & time == "144h") | (genotype == "wt" & rep %in% c("rep3", "rep4") & time != "144h")))
to.plot.2$gene <- gsub("Hox", "", gsub("_ratio", "", to.plot.2$gene))
to.plot.2$gene <- factor(to.plot.2$gene, levels = c("d13", "d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.2$time <- factor(to.plot.2$time, levels = c("96h", "120h", "144h"))

# Manual filtering:
to.plot.2$ratio_wt[to.plot.2$time != "144h" & to.plot.2$gene == "d13"] <- NA
to.plot.2$ratio_wt[to.plot.2$time == "96h" & to.plot.2$gene %in% c("d10", "d11")] <- NA

figures <- list("Fig7B" = list(geno = "Del(CBS1-5)",
                                scales = list(scale_y_continuous(breaks = c(0, 1, 2, 4, 6),
                                                                 limits = c(0, 6),
                                                                 expand = c(0, 0)),
                                              scale_y_continuous(breaks = c(5, 10, 15, 20),
                                                                 limits = c(4, 20),
                                                                 expand = c(0, 0)),
                                              scale_y_continuous(breaks = c(0, 1, 2, 4),
                                                                 limits = c(0, 4),
                                                                 expand = c(0, 0)),
                                              scale_y_continuous(breaks = c(5, 10, 15),
                                                                 limits = c(5, 15),
                                                                 expand = c(0, 0)),
                                              scale_y_continuous(breaks = c(0, 1, 2, 4),
                                                                 limits = c(0, 5),
                                                                 expand = c(0, 0)))))

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.2, genotype %in% c("wt", figures[[fig.name]][["geno"]]) & is.finite(ratio_wt))
  if (fig.name == "Fig7B") {
    df.stats$time <- as.character(df.stats$time)
    temp.df <- subset(df.stats, time != "96h")
    df.stats$time[df.stats$time != "96h"] <- paste("high", df.stats$time[df.stats$time != "96h"])
    temp.df$time <- paste("low", temp.df$time)
    df.stats <- rbind(df.stats, temp.df)
    df.stats$time <- factor(df.stats$time, levels = c("96h", "high 120h", "low 120h", "high 144h", "low 144h"))
  }
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  ns <- table(df.stats$split)
  df.stats.split <- split(subset(df.stats, split %in% names(ns[ns > 3])), factor(subset(df.stats, split %in% names(ns[ns > 3]))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("Normalized FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, ".pdf")), g, height = 5, width = 7)
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}

# No ratio
to.plot.2.noratio <- subset(all.FPKM.values.v2, 
                            gene %in% paste0("Hox", c("d13", "d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1")) &
                              (genotype == "Del(CBS1-5)"  | (genotype == "wt" & rep %in% c("rep1", "rep2") & time == "144h") | (genotype == "wt" & rep %in% c("rep3", "rep4") & time != "144h")))
to.plot.2.noratio$gene <- gsub("Hox", "", to.plot.2.noratio$gene)
to.plot.2.noratio$gene <- factor(to.plot.2.noratio$gene, levels = c("d13", "d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.2.noratio$time <- factor(to.plot.2.noratio$time, levels = c("96h", "120h", "144h"))

# Manual filtering:
to.plot.2.noratio$ratio_wt[to.plot.2.noratio$time != "144h" & to.plot.2.noratio$gene == "d13"] <- NA
to.plot.2.noratio$ratio_wt[to.plot.2.noratio$time == "96h" & to.plot.2.noratio$gene %in% c("d10", "d11")] <- NA

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.2.noratio, genotype %in% c("wt", figures[[fig.name]][["geno"]]) & is.finite(ratio_wt))
  if (fig.name == "Fig7B") {
    df.stats$time <- as.character(df.stats$time)
    temp.df <- subset(df.stats, time != "96h")
    df.stats$time[df.stats$time != "96h"] <- paste("high", df.stats$time[df.stats$time != "96h"])
    temp.df$time <- paste("low", temp.df$time)
    df.stats <- rbind(df.stats, temp.df)
    df.stats$time <- factor(df.stats$time, levels = c("96h", "high 120h", "low 120h", "high 144h", "low 144h"))
  }
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  ns <- table(df.stats$split)
  df.stats.split <- split(subset(df.stats, split %in% names(ns[ns > 3])), factor(subset(df.stats, split %in% names(ns[ns > 3]))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, "_not_norm.pdf")), g, height = 5, width = 7)
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "_not_norm_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}

to.plot.3 <- subset(all.FPKM.values.v2, 
                    gene %in% paste0("Hox", c("d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"), "_ratio") &
                      genotype %in% c("Del(CBS4)", "wt") & rep %in% c("rep1", "rep2") & time != "96h")
to.plot.3$gene <- gsub("Hox", "", gsub("_ratio", "", to.plot.3$gene))
to.plot.3$gene <- factor(to.plot.3$gene, levels = c("d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.3$time <- factor(to.plot.3$time, levels = c("96h", "120h", "144h"))

# Manual filtering:
to.plot.3$ratio_wt[to.plot.3$time == "120h" & to.plot.3$gene %in% c("d10", "d11")] <- NA

figures <- list("FigS10B" = list(geno = "Del(CBS4)",
                               scales = list(scale_y_continuous(breaks = c(0, 1, 2),
                                                                limits = c(0, 2),
                                                                expand = c(0, 0)),
                                             scale_y_continuous(breaks = 0:4,
                                                                limits = c(0, 4),
                                                                expand = c(0, 0)))))

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.3, genotype %in% c("wt", figures[[fig.name]][["geno"]]))
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  df.stats.split <- split(subset(df.stats, time != "120h" | !gene %in% c("d10", "d11")), factor(subset(df.stats, time != "120h" | !gene %in% c("d10", "d11"))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("Normalized FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, ".pdf")), g, height = 5, width = 7)
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}

to.plot.3.noratio <- subset(all.FPKM.values.v2, 
                            gene %in% paste0("Hox", c("d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1")) &
                              genotype %in% c("Del(CBS4)", "wt") & rep %in% c("rep1", "rep2") & time != "96h")
to.plot.3.noratio$gene <- gsub("Hox", "", to.plot.3.noratio$gene)
to.plot.3.noratio$gene <- factor(to.plot.3.noratio$gene, levels = c("d11", "d10", "d9", "d8", "d4", "d3", "d1", "a9", "a5", "a1"))
to.plot.3.noratio$time <- factor(to.plot.3.noratio$time, levels = c("96h", "120h", "144h"))

# Manual filtering:
to.plot.3.noratio$ratio_wt[to.plot.3.noratio$time == "120h" & to.plot.3.noratio$gene %in% c("d10", "d11")] <- NA

for (fig.name in names(figures)) {
  df.stats <- subset(to.plot.3.noratio, genotype %in% c("wt", figures[[fig.name]][["geno"]]))
  df.stats$split <- paste0(df.stats$time, "_", df.stats$gene)
  df.stats.split <- split(subset(df.stats, time != "120h" | !gene %in% c("d10", "d11")), factor(subset(df.stats, time != "120h" | !gene %in% c("d10", "d11"))$split))
  tests <- lapply(df.stats.split, function(df){t.test(ratio_wt ~ genotype, df, paired = FALSE, var.equal = FALSE)})
  p.values <- sapply(tests, "[[", "p.value")
  temp.df <- unique(df.stats[, c("time", "gene", "split")])
  temp.df$pvals <- p.values[temp.df$split]
  temp.df$label <- "NS"
  temp.df$label[temp.df$pvals < 0.05] <- "*"
  temp.df$label[temp.df$pvals < 0.01] <- "**"
  temp.df$label[temp.df$pvals < 0.001] <- "***"
  temp.df$label[temp.df$pvals < 0.0001] <- "****"
  temp.df$label[is.na(temp.df$pvals)] <- "NA"
  working.df <- subset(df.stats, genotype != "wt")
  max.val.df <- ddply(working.df , .(gene, time), summarize, ratio_wt = 1.1 * max(ratio_wt))
  temp.df <- merge(temp.df, max.val.df)
  g <- ggplot(working.df, aes(x = gene, y = ratio_wt, color = gene, shape = gene)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "bar", fill = NA, width = 0.6) + 
    stat_summary(fun.data = simple.fun, 
                 geom = "errorbar", width = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_text(data = temp.df, aes(label = label), color = "black") +
    facet_grid(time ~ ., scale="free_y") +
    theme_classic() +
    scale_shape_manual(values = gene.shapes) +
    scale_color_manual(values = gene.colors) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0) +
    ylab("FPKM vs wt") +
    ggtitle(figures[[fig.name]][["geno"]]) +
    facetted_pos_scales(y = figures[[fig.name]][["scales"]])
  ggsave(file.path(directory.with.plots, paste0(fig.name, "_not_norm.pdf")), g, height = 5, width = 7)
  
  
  # Export data:
  write.table(cast(df.stats[, c("rep", "time", "gene", "genotype", "ratio_wt")],
                   rep ~ time + gene + genotype, value = "ratio_wt"),
              file.path(directory.with.data, paste0(fig.name, "_not_norm_plotted_data.txt")),
              sep = "\t", quote = F, row.names = F)
}
