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
if (!"FSA" %in% installed.packages()) {
  install.packages("FSA", repos = "https://stat.ethz.ch/CRAN/")
}
library(FSA)
if (!"ggpubr" %in% installed.packages()) {
  install.packages("ggpubr", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggpubr)
if (!"rstatix" %in% installed.packages()) {
  install.packages("rstatix", repos = "https://stat.ethz.ch/CRAN/")
}
library(rstatix)

# Function to plot error bars
min.max <- function(v){return(list(ymin = min(v),
                                   ymax = max(v)))}

# Set paths
directory.with.v4C <- commandArgs(TRUE)[1]
directory.with.github <- commandArgs(TRUE)[2]
directory.with.plots <- file.path(directory.with.github, "cHi-C", "plots")
directory.with.data <- file.path(directory.with.plots, "..", "raw_values")
# Create non-existing output directories
if (!dir.exists(directory.with.plots)) {
  dir.create(directory.with.plots, recursive = T, showWarnings = F)
}
if (!dir.exists(directory.with.data)) {
  dir.create(directory.with.data, recursive = T, showWarnings = F)
}

regions.df <- read.delim(file.path(directory.with.github, "annotations", "v4C_quantified.bed"), header = F)
regions <- apply(regions.df, 1, function(v){seq(as.integer(v[2]), as.integer(v[3]) - 1, 1000)})
names(regions) <- regions.df[, 4]
# Manual:
regions[["sub-TAD1-noCTCF"]] <- c(regions[["5p-sub-TAD1"]], regions[["middle-sub-TAD1"]], regions[["3p-sub-TAD1"]])
regions[["CBS-sub-TAD1"]] <- c(regions[["CBS-sub-TAD1-1"]], regions[["CBS-sub-TAD1-2"]], regions[["CBS-CS38-40"]])

figures <- list("FigS6A" = list("vp" = paste0("CBS", 6:1),
                                "time" = 48,
                                "genotype" = "wt",
                                "x" = "vp",
                                "color" = "time",
                                "regions" = "sub-TAD1-noCTCF",
                                "title" = "sub-TAD1-noCTCF"),
                "FigS6B" = list("vp" = paste0("CBS", 6:1),
                                "time" = 48,
                                "genotype" = "wt",
                                "x" = "vp",
                                "color" = "time",
                                "regions" = "CBS-sub-TAD1",
                                "title" = "CBS-sub-TAD1"),
                "FigS6C-1" = list("vp" = paste0("CBS", 6:1),
                                  "time" = 48,
                                  "genotype" = "wt",
                                  "x" = "vp",
                                  "color" = "time",
                                  "regions" = "CBS-CS38-40",
                                  "title" = "CBS-CS38-40"),
                "FigS6C-2" = list("vp" = paste0("CBS", 6:1),
                                  "time" = 96,
                                  "genotype" = "wt",
                                  "x" = "vp",
                                  "color" = "time",
                                  "regions" = "CBS-CS38-40",
                                  "title" = "CBS-CS38-40"),
                "FigS6D" = list("vp" = "d1-d4",
                                "time" = seq(48, 168, 24),
                                "genotype" = "wt",
                                "x" = "time",
                                "color" = "time",
                                "regions" = "sub-TAD1-noCTCF",
                                "title" = "d1-d4 vs sub-TAD1-noCTCF"),
                "FigS6E-1" = list("vp" = "d9",
                                  "time" = seq(48, 168, 24),
                                  "genotype" = "wt",
                                  "x" = "time",
                                  "color" = "time",
                                  "regions" = "sub-TAD1-noCTCF",
                                  "title" = "d9 vs subTAD1-noCTCF"),
                "FigS6E-2" = list("vp" = "d9",
                                  "time" = seq(48, 168, 24),
                                  "genotype" = "wt",
                                  "x" = "time",
                                  "color" = "time",
                                  "regions" = "3p-sub-TAD1",
                                  "title" = "d9 vs 3p sub-TAD1"),
                "FigS8A" = list("vp" = "d11",
                                "time" = seq(48, 168, 48),
                                "genotype" = "wt",
                                "x" = "time",
                                "color" = "time",
                                "regions" = "Evx2-Hoxd12",
                                "test" = "kd",
                                "group" = "time",
                                "ref" = "48",
                                "title" = "d11 vs Evx2-d12"),
                "FigS8B" = list("vp" = "d9",
                                "time" = seq(48, 168, 48),
                                "genotype" = "wt",
                                "x" = "time",
                                "color" = "time",
                                "regions" = "Evx2-Hoxd12",
                                "test" = "kd",
                                "group" = "time",
                                "ref" = "48",
                                "title" = "d9 vs Evx2-d12"),
                "FigS8C" = list("vp" = "d8",
                                "time" = seq(48, 168, 48),
                                "genotype" = "wt",
                                "x" = "time",
                                "color" = "time",
                                "regions" = "Evx2-Hoxd12",
                                "test" = "kd",
                                "group" = "time",
                                "ref" = "48",
                                "title" = "d8 vs Evx2-d12"),
                "FigS14A-1" = list("vp" = paste0("CBS", c(5, 4, 2, 1)),
                                  "time" = 96,
                                  "genotype" = "wt",
                                  "x" = "vp",
                                  "color" = "genotype",
                                  "regions" = "CBS-CS38-40",
                                  "test" = "kd",
                                  "group" = "vp",
                                  "title" = "CBS-CS38-40"),
                "FigS14A-2" = list("vp" = paste0("CBS", c(5, 4, 2, 1)),
                                  "time" = 96,
                                  "genotype" = "Del(d1-d4)",
                                  "x" = "vp",
                                  "color" = "genotype",
                                  "regions" = "CBS-CS38-40",
                                  "test" = "kd",
                                  "group" = "vp",
                                  "title" = "CBS-CS38-40"),
                "FigS14B" = list("vp" = paste0("d", c(13, 11, 9, 8)),
                                "time" = 96,
                                "genotype" = c("wt", "Del(sub-TAD1)"),
                                "x" = "vp",
                                "color" = "genotype",
                                "regions" = "CBS-CS38-40",
                                "test" = "wilcox",
                                "group" = "genotype",
                                "split" = "vp",
                                "ref" = "wt",
                                "title" = "CBS-CS38-40"),
                "FigS15B" = list("vp" = paste0("d", c(13, 11, 9)),
                                 "time" = 96,
                                 "genotype" = c("wt", "Del(CBS1-5)"),
                                 "x" = "vp",
                                 "color" = "genotype",
                                 "regions" = "Hoxd3-4",
                                 "test" = "wilcox",
                                 "group" = "genotype",
                                 "split" = "vp",
                                 "ref" = "wt",
                                 "title" = "CBS-CS38-40"),
                "FigS16C" = list("vp" = paste0("d", c(13, 11, 9)),
                                 "time" = 96,
                                 "genotype" = c("wt", "Del(CBS1-5)"),
                                 "x" = "vp",
                                 "color" = "genotype",
                                 "regions" = "5p-sub-TAD1",
                                 "test" = "wilcox",
                                 "group" = "genotype",
                                 "split" = "vp",
                                 "ref" = "wt",
                                 "title" = "5p sub-TAD1")
)

fixed.colors <- c("48" = "#3d58a7", "72" = "#ee3124", "96" = "#4ab848", "120" = "#9753a1", "144" = "#f7931d", "168" = "#040707",
                  "wt" = "#4f9041", "Del(d1-d4)" = "#ee2f24", "Del(sub-TAD1)" = "#ee2f24", "Del(CBS1-5)" = "#ee2f24"
)

all.bedgraphs.fn <- list.files(path = directory.with.v4C, pattern = ".bedgraph", full.names = T)
all.bedgraphs.df <- lapply(all.bedgraphs.fn, function(fn) {
  df <- read.delim(fn, header = F)
  colnames(df) <- c("chr", "start", "end", "value")
  df$sample <- gsub(".bedgraph$", "", basename(fn))
  return(df)
})
all.bedgraphs.values <- do.call(rbind, all.bedgraphs.df)


samplesplan <- data.frame(sample = unique(as.character(all.bedgraphs.values$sample)))
samplesplan$genotype <- sapply(strsplit(samplesplan$sample, "_"), "[[", 1)
samplesplan$time <- as.numeric(gsub("h$", "", sapply(strsplit(samplesplan$sample, "_"), "[[", 2)))
samplesplan$vp <- sapply(strsplit(samplesplan$sample, "_"), tail, n = 1)

input.data <- merge(all.bedgraphs.values, samplesplan)

for (fig in names(figures)) {
  working.df <- subset(input.data, vp %in% figures[[fig]][["vp"]] &
                         time %in% figures[[fig]][["time"]] &
                         genotype %in% figures[[fig]][["genotype"]] &
                         start %in% regions[[figures[[fig]][["regions"]]]])
  ylab <- "Virtual4C contacts"
  if ("Del(d1-d4)" %in% working.df$genotype) {
    # The region need to be shifted
    working.df <- subset(working.df, genotype != "Del(d1-d4)")
    my.starts <- regions[[figures[[fig]][["regions"]]]]
    my.shifted.starts <- c(my.starts[my.starts < 74719771], my.starts[my.starts > 74765994] - 45000)
    working.df <- rbind(working.df,
                        subset(input.data, vp %in% figures[[fig]][["vp"]] &
                                 time %in% figures[[fig]][["time"]] &
                                 genotype == "Del(d1-d4)" &
                                 start %in% my.shifted.starts))
  }
  if ("Del(sub-TAD1)" %in% working.df$genotype) {
    # The region need to be shifted
    working.df <- subset(working.df, genotype != "Del(sub-TAD1)")
    my.starts <- regions[[figures[[fig]][["regions"]]]]
    my.shifted.starts <- c(my.starts[my.starts < 74768588], my.starts[my.starts > 75133800] - 365000)
    working.df <- rbind(working.df,
                        subset(input.data, vp %in% figures[[fig]][["vp"]] &
                                 time %in% figures[[fig]][["time"]] &
                                 genotype == "Del(sub-TAD1)" &
                                 start %in% my.shifted.starts))
  }
  if (grepl("FigS8", fig)) {
    # Normalize to average of 48h
    working.df$original.value <- working.df$value
    working.df$value <- working.df$value / mean(working.df$value[working.df$time == 48])
    ylab <- "Normalized Virtual4C contacts"
  }
  working.df$time <- factor(working.df$time, levels = figures[[fig]][["time"]])
  working.df$vp <- factor(working.df$vp, levels = figures[[fig]][["vp"]])
  working.df$genotype <- factor(working.df$genotype, levels = figures[[fig]][["genotype"]])
  sub.title <- ""
  results <- data.frame(label = character(0), value = numeric(0))
  results[, figures[[fig]][["x"]]] <- character(0)
  results[, figures[[fig]][["color"]]] <- character(0)
  stat.test <- NULL
  if ("test" %in% names(figures[[fig]]) && figures[[fig]][["test"]] == "kd") {
    if ("split" %in% names(figures[[fig]])) {
      my.working.df <- split(working.df, factor(working.df[, figures[[fig]][["split"]]]))
    } else {
      my.working.df <- list('all' = working.df)
    }
    for (my.group in names(my.working.df)) {
      p.val <- kruskal.test(x = my.working.df[[my.group]]$value, g = my.working.df[[my.group]][, figures[[fig]][["group"]]])$p.value
      if (p.val < 0.05) {
        sub.title <- paste0(sub.title, "\n", my.group, ": Kruskal: p-value = ", signif(p.val, 2))
        # For the stars with ref:
        if ("ref" %in% names(figures[[fig]])) {
          temp.results <- dunnTest(x = my.working.df[[my.group]]$value,
                                   g = my.working.df[[my.group]][, figures[[fig]][["group"]]],
                                   method = "none")$res
          # Only keep comparisons involving ref
          temp.results <- subset(temp.results, grepl(paste0("^", figures[[fig]][["ref"]], " - "), Comparison) |
                                   grepl(paste0(" - ", figures[[fig]][["ref"]], "$"), Comparison))
          # Re-evaluate padj:
          temp.results$P.adj <- p.adjust(temp.results$P.unadj, method = "bonferroni")
          temp.results$label <- "NS"
          temp.results$label[temp.results$P.adj < 0.05] <- "*"
          temp.results$label[temp.results$P.adj < 0.01] <- "**"
          temp.results$label[temp.results$P.adj < 0.001] <- "***"
          temp.results$label[temp.results$P.adj < 0.0001] <- "****"
          temp.results[, figures[[fig]][["group"]]] <- gsub(paste0("^", figures[[fig]][["ref"]], " - "), "", temp.results$Comparison)
          temp.results[, figures[[fig]][["group"]]] <- gsub(paste0(" - ", figures[[fig]][["ref"]], "$"), "", temp.results[, figures[[fig]][["group"]]])
          max.val <- ddply(my.working.df[[my.group]], figures[[fig]][["group"]], summarise, max.val = max(value))
          temp.results$value <- 1.05 * max.val$max.val[match(temp.results[, figures[[fig]][["group"]]], max.val[, figures[[fig]][["group"]]])]
          if (!all(colnames(results) %in% colnames(temp.results))) {
            for (cn in setdiff(colnames(results), colnames(temp.results))) {
              temp.results[, cn] <- unique(my.working.df[[my.group]][, cn])
            }
          }
          results <- rbind(results, temp.results[, colnames(results)])
        }
        # For the brackets:
        stat.test <- dunn_test(data = my.working.df[[my.group]],
                               formula = as.formula(paste("value ~", figures[[fig]][["group"]])),
                               p.adjust.method = "bonferroni")
        # Remove the ns:
        stat.test <- subset(stat.test, p.adj.signif != "ns")
        stat.test <- stat.test %>%
          add_y_position()
      } else {
        sub.title <- paste0(sub.title, "\n", my.group, ": Kruskal NS")
      }
    }
    if (grepl("p-value", sub.title) & "ref" %in% names(figures[[fig]])) {
      sub.title <- paste0(sub.title, "\ncolored asterisks corresponds to Dunn's adj only vs ", figures[[fig]][["ref"]])
    }
  }
  if ("test" %in% names(figures[[fig]]) && figures[[fig]][["test"]] == "wilcox") {
    if ("split" %in% names(figures[[fig]])) {
      my.working.df <- split(working.df, factor(working.df[, figures[[fig]][["split"]]]))
    } else {
      my.working.df <- list('all' = working.df)
    }
    pvals <- lapply(my.working.df, function(df){wilcox.test(subset(df, df[, figures[[fig]][["group"]]] == figures[[fig]][["ref"]])$value,
                                                            subset(df, df[, figures[[fig]][["group"]]] != figures[[fig]][["ref"]])$value,
                                                            paired = FALSE)$p.value})
    max.val <- sapply(lapply(my.working.df, "[[", "value"), max)
    temp.results <- data.frame(pvals = unlist(pvals), value = 1.05 * unlist(max.val))
    if ("split" %in% names(figures[[fig]])) {
      temp.results[, figures[[fig]][["split"]]] <- names(pvals)
    }
    non.ref <- lapply(my.working.df, function(df){
      unique(as.character(subset(df, df[, figures[[fig]][["group"]]] != figures[[fig]][["ref"]])[, figures[[fig]][["group"]]]))
    })
    temp2.res <- data.frame(group = unlist(non.ref), split = rep(names(non.ref), sapply(non.ref, length)))
    colnames(temp2.res) <- c(figures[[fig]][["group"]], figures[[fig]][["split"]])
    temp.results <- merge(temp.results, temp2.res)
    temp.results$label <- "NS"
    temp.results$label[temp.results$pvals < 0.05] <- "*"
    temp.results$label[temp.results$pvals < 0.01] <- "**"
    temp.results$label[temp.results$pvals < 0.001] <- "***"
    temp.results$label[temp.results$pvals < 0.0001] <- "****"
    if (!all(colnames(results) %in% colnames(temp.results))) {
      for (cn in setdiff(colnames(results), colnames(temp.results))) {
        temp.results[, cn] <- unique(my.working.df[[my.group]][, cn])
      }
    }
    results <- temp.results[, colnames(results)]
    sub.title <- paste0("asterisks corresponds to Mann-Whitney vs ", figures[[fig]][["ref"]])
  }
  # For reproducibility
  set.seed(1)
  max.value <- max(working.df$value)
  if (!is.null(stat.test)) {
    max.value <- max(max.value, max(stat.test$y.position))
  }
  # I use as error bars the min max
  g <- ggplot(working.df, aes_string(x = figures[[fig]][["x"]], y = "value", color = figures[[fig]][["color"]])) +
    stat_summary(fun.data = min.max, 
                 geom = "errorbar",
                 position = position_dodge(0.8)) +
    geom_boxplot(fill = "white",
                 position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(0.8)) +
    geom_text(data = results, aes(label = label)) +
    theme_classic() +
    ylab(ylab) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max.value * 1.1)) +
    scale_color_manual(values = fixed.colors[as.character(unique(working.df[, figures[[fig]][["color"]]]))]) +
    labs(title = figures[[fig]][["title"]],
         caption = sub.title)
  
  if (!is.null(stat.test)) {
    g <- g +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
  }
  ggsave(file.path(directory.with.plots, paste0(fig, ".pdf")), height = 5, width = 5)
  # Export data:
  write.table(cast(working.df, chr + start + end ~ sample), file.path(directory.with.data, paste0(fig, "_raw_data.txt")), sep = "\t", quote = F, row.names = F)
}

