# Load packages
if (!"reshape" %in% installed.packages()) {
  install.packages("reshape", repos = "https://stat.ethz.ch/CRAN/")
}
library(reshape)
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
if (!"plyr" %in% installed.packages()) {
  install.packages("plyr", repos = "https://stat.ethz.ch/CRAN/")
}
library(plyr)

# Set paths
directory.with.quantifs <- commandArgs(TRUE)[1]
directory.with.plots <- commandArgs(TRUE)[2]
# Create non-existing output directory
if (!dir.exists(directory.with.plots)) {
  dir.create(directory.with.plots, recursive = T, showWarnings = F)
}

figures.type <- list('Fig1D' = 'coo',
                     'Fig2B' = 'coo',
                     'FigS1' = 'time',
                     'FigS4A' = 'coo',
                     'FigS4B' = 'time',
                     'FigS5B' = 'time')
figures.group <- list('FigS1' = list('heatmap' = list("1" = "HoxD",
                                                      "2" = c("C-DOM", "SubTAD2", "SubTAD1","T-DOM")),
                                     'points' = list("3" = rev(c("T6: CS5", "T5: CS85-31", "T4: CS31-32", "T3: CS34", "T2: CS37", "T1: CS39")),
                                                     "4" = rev(c("C4: IslandI", "C3: IslandII-III", "C2: IslandIII", "C1: GT2")))
),
'FigS4B' = list('points' = list("1" = paste0("CBS", c(1, 2, 4, 5)),
                                "2" = paste0("CBS", c(3, 6:9)))
),
'FigS5B' = list('points' = list("1" = paste0("CD-CBS", 1:5),
                                "2" = paste0("TD-CBS", 1:9))
)
)
ratios <- data.frame(protein = character(0), time = character(0), ratio = numeric(0))
for (fig in names(figures.type)) {
  df <- read.delim(file.path(directory.with.quantifs, paste0(fig, ".txt")))
  
  # I process the colnames of the data frame:
  colnames(df) <- gsub("^X\\.", "", colnames(df))
  colnames(df) <- gsub("^\\.", "", colnames(df))
  colnames(df) <- gsub("\\.$", "", colnames(df))
  colnames(df) <- gsub("_Normalized.bigwig$", "", colnames(df))
  
  if (nrow(df) == 10) {
    df$region <- factor(paste0("R_", 10:1), levels = paste0("R_", 10:1))
  } else {
    region.df <- read.delim(file.path(directory.with.quantifs, paste0('regions_', gsub("Fig", "", fig), '.bed')), header = F)[, 1:4]
    colnames(region.df) <- c("chr", "start", "end", "region")
    region.df$region <- factor(region.df$region, levels = region.df$region)
    df <- merge(df, region.df)
  }
  all.values <- melt(df, id.vars = c("chr", "start", "end", "region"), variable_name = "sample")
  
  samplesplan <- data.frame(sample = unique(as.character(all.values$sample)))
  samplesplan$time <- sapply(strsplit(samplesplan$sample, "_"), "[[", 2)
  
  protein <- strsplit(samplesplan$sample, "_")[[1]][3]
  # Put space between l and I:
  protein <- gsub("lI", "l I", protein)
  
  input.data <- merge(all.values, samplesplan)
  input.data$time <- factor(input.data$time, levels = paste0(seq(168, 48, -12), "h"))
  input.data$time.num <- as.numeric(gsub("h", "", as.character(input.data$time)))
  all.times <- data.frame(time.num = sort(unique(input.data$time.num)))
  all.times.inter <- diff(all.times$time.num)
  all.times$before <- NA
  all.times$after <- NA
  all.times$before[2:nrow(all.times)] <- all.times$time.num[2:nrow(all.times)] - all.times.inter / 2
  all.times$after[1:(nrow(all.times) - 1)] <- all.times$time.num[1:(nrow(all.times) - 1)] + all.times.inter / 2
  all.times$before[1] <- all.times$time.num[1] - median(all.times.inter) / 2
  all.times$after[nrow(all.times)] <- all.times$time.num[nrow(all.times)] + median(all.times.inter) / 2
  input.data <- merge(input.data, all.times)
  
  if (fig == "FigS5B") {
    # Normalize to 72h
    input.data <- ddply(input.data, .(chr, start, region), mutate,
                        new.value = value / value[time == "72h"])
    input.data$original.value <- input.data$value
    input.data$value <- input.data$new.value
  }
  if (figures.type[fig] == "coo") {
    if (fig == "FigS4A") {
      ggplot(subset(input.data, time != "48h"), aes(x = region, y = time, fill = value)) +
        geom_tile() +
        scale_fill_distiller(paste("Coverage of", protein),
                             palette = "OrRd",
                             direction = 0) +
        xlab("") +
        ggtitle(protein) +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_blank(),
              strip.background = element_rect(fill = NA, colour = NA))
      ggsave(file.path(directory.with.plots, paste0(fig, ".pdf")), height = 3, width = 5)
    } else {
      highlighted.time.point <- 96
      before <- all.times$before[all.times$time.num == highlighted.time.point]
      after <- all.times$after[all.times$time.num == highlighted.time.point]
      highlight <- data.frame(time = c(before, before, after, after, before),
                              x = c(74707561, 74697514, 74697514, 74707561, 74707561))
      ggplot(input.data) +
        geom_rect(aes(xmin = start, xmax = end, ymin = before, ymax = after, fill = value)) +
        scale_fill_distiller("Coverage",
                             palette = "OrRd",
                             direction = 0) +
        geom_path(data = highlight, aes(x, time), lty = 2) +
        xlab("") +
        ggtitle(protein) +
        theme_classic() +
        scale_y_reverse("time",
                        breaks = unique(input.data$time.num),
                        labels = unique(as.character(input.data$time)),
                        expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))
      ggsave(file.path(directory.with.plots, paste0(fig, ".pdf")), height = 3, width = 5)
    }
    if ("R_7" %in% input.data$region & protein %in% c("H3K27ac", "Pol II")) {
      ratios <- rbind(ratios,
                      data.frame(protein = protein, time = c("96h", "120h"),
                                 ratio = c(input.data$value[input.data$time == "96h" & input.data$region == "R_7"] / 
                                             input.data$value[input.data$time == "96h" & input.data$region == "R_4"],
                                           input.data$value[input.data$time == "120h" & input.data$region == "R_7"] / 
                                             input.data$value[input.data$time == "120h" & input.data$region == "R_4"])))
    }
  } else if (figures.type[fig] == "time") {
    # First the colormaps
    for (i in names(figures.group[[fig]][["heatmap"]])) {
      my.regions <- figures.group[[fig]][["heatmap"]][[i]]
      working.df <- subset(input.data, region %in% my.regions & time != "48h")
      all.times <- data.frame(time.num = sort(unique(working.df$time.num)))
      all.times.inter <- diff(all.times$time.num)
      all.times$before <- NA
      all.times$after <- NA
      all.times$before[2:nrow(all.times)] <- all.times$time.num[2:nrow(all.times)] - all.times.inter / 2
      all.times$after[1:(nrow(all.times) - 1)] <- all.times$time.num[1:(nrow(all.times) - 1)] + all.times.inter / 2
      all.times$before[1] <- all.times$time.num[1] - median(all.times.inter) / 2
      all.times$after[nrow(all.times)] <- all.times$time.num[nrow(all.times)] + median(all.times.inter) / 2
      working.df <- merge(working.df, all.times)
      y.center <- seq(1, length(my.regions))
      names(y.center) <- my.regions
      working.df$ymin <- y.center[as.character(working.df$region)] - 0.5
      working.df$ymax <- y.center[as.character(working.df$region)] + 0.5
      ggplot(working.df, aes(xmin = before, xmax = after, ymin = ymin, ymax = ymax, fill = value)) +
        geom_rect() +
        scale_fill_distiller("Coverage",
                             palette = "OrRd",
                             direction = 0) +
        xlab("") +
        theme_classic() +
        scale_x_continuous("time",
                           breaks = unique(working.df$time.num),
                           labels = unique(as.character(working.df$time)),
                           expand = c(0, 0)) +
        scale_y_continuous("", breaks = y.center, labels = names(y.center),
                           expand = c(0, 0))
      ggsave(file.path(directory.with.plots, paste0(fig, "-", i, ".pdf")), height = 0.5 * (1 + length(my.regions)), width = 5)
    }
    # Then the lines
    for (i in names(figures.group[[fig]][["points"]])) {
      my.regions <- figures.group[[fig]][["points"]][[i]]
      working.df <- subset(input.data, region %in% my.regions & time != "48h")
      working.df$region <- factor(working.df$region, levels = my.regions)
      title <- paste("Coverage of", protein)
      if (fig == "FigS5B"){
        title <- paste(title, "Normalized to 72h")
      }
      ggplot(working.df, aes(x = time.num, y = value, color = region, shape = region)) +
        geom_point() +
        stat_summary(fun = mean, aes(group = region), 
                     geom = "line") +
        ylab(title) +
        theme_classic() +
        scale_x_continuous("time",
                           breaks = unique(working.df$time.num),
                           labels = unique(as.character(working.df$time)))
      ggsave(file.path(directory.with.plots, paste0(fig, "-", i, ".pdf")), height = 3, width = 5)
    }
  }
}
# Plot the ratios in Fig2C left
ratios$time <- factor(ratios$time, levels = c("96h", "120h"))
ggplot(ratios, aes(x = time, y = ratio)) +
  geom_line(aes(group = protein)) +
  facet_grid(. ~ protein) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("") +
  ylab("d9/d4 coverage ratio") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))
ggsave(file.path(directory.with.plots, paste0("Fig2Ca.pdf")), height = 3, width = 3)
