# @author: Alexandre Mayran
# Linting by Lucille Delisle

# This analysis was initially run by Alexandre with:
# > sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=French_Switzerland.1252  LC_CTYPE=French_Switzerland.1252    LC_MONETARY=French_Switzerland.1252
# [4] LC_NUMERIC=C                        LC_TIME=French_Switzerland.1252
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
# [1] BiocManager_1.30.10         ggpubr_0.4.0                tidyr_1.1.2
# [4] RColorBrewer_1.1-2          scales_1.1.1                MAST_1.16.0
# [7] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0
# [10] GenomicRanges_1.42.0        GenomeInfoDb_1.26.2         IRanges_2.24.1
# [13] S4Vectors_0.28.1            BiocGenerics_0.36.0         MatrixGenerics_1.2.1
# [16] matrixStats_0.58.0          sctransform_0.3.2           ggplot2_3.3.5
# [19] SeuratObject_4.0.0          Seurat_4.0.2                dplyr_1.0.3
# [22] usefulLDfunctions_0.1.5

# And confirmed by Lucille with:
# > sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Arch Linux
# 
# Matrix products: default
# BLAS:   /usr/lib/libblas.so.3.10.1
# LAPACK: /usr/lib/liblapack.so.3.10.1
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=C                  LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
# [1] RColorBrewer_1.1-3          scales_1.2.0                MAST_1.22.0                 SingleCellExperiment_1.18.0
# [5] SummarizedExperiment_1.26.0 Biobase_2.56.0              GenomicRanges_1.48.0        GenomeInfoDb_1.32.1
# [9] IRanges_2.30.0              S4Vectors_0.34.0            BiocGenerics_0.42.0         MatrixGenerics_1.8.0
# [13] matrixStats_0.62.0          sctransform_0.3.3           ggplot2_3.3.5               SeuratObject_4.0.4
# [17] Seurat_4.1.0                dplyr_1.0.9                 usefulLDfunctions_0.1.6

# Also run a third time on another machine, this is what has been used to generate what is on GEO in on GitHub:
# > sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.3 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2          scales_1.2.1                MAST_1.18.0                 SingleCellExperiment_1.14.1
# [5] SummarizedExperiment_1.22.0 MatrixGenerics_1.4.3        matrixStats_0.61.0          sctransform_0.3.3          
# [9] dplyr_1.0.9                 MuSiC_0.2.0                 nnls_1.4                    ggh4x_0.2.3                
# [13] plyr_1.8.6                  DWLS_0.1.0                  ggpubr_0.4.0                ggplot2_3.4.0              
# [17] reshape2_1.4.4              rtracklayer_1.52.1          GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
# [21] IRanges_2.26.0              S4Vectors_0.30.2            Biobase_2.52.0              BiocGenerics_0.38.0        
# [25] sp_1.4-7                    SeuratObject_4.1.0          Seurat_4.1.1                usefulLDfunctions_0.1.6    

############
# These 2 lines are specific to windows:
# Define library directory
.libPaths(c("library", "C:/Users/Mayran/Desktop/lab/R/2020.analyses/AM.analyses/Program FilesRR-3.6.1/library/"))

#to avoid reaching memory limits define higher memory limits
memory.size(50000)
############

# Install required packages
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
if (!"usefulLDfunctions" %in% installed.packages()) {
  devtools::install_github("lldelisle/usefulLDfunctions")
}
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("sctransform")
safelyLoadAPackageInCRANorBioconductor("MAST")
safelyLoadAPackageInCRANorBioconductor("scales")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")

# Introduce variable that will be used
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
output.directory <- "WT/output.files/"
output.metrics <- paste0(output.directory, "metrics/")
output.matrices <- paste0(output.directory, "matrices/")

# set working directory
setwd()

# Fill with the github path
github.directory <- ""
# Fill with path with matrices
input.directory <- ""

# create lists that will be filled
# all.sorted.raw.data contains the output of Read10x
all.sorted.raw.data <- list()
# seurat.all.sorted contains the output of CreateSeuratObject with minimal filtering
seurat.all.sorted <- list()
# seurat.all.sorted.subset contains cells filtered for low quality cells/doublets
seurat.all.sorted.subset <- list()
# mean.rna contains the average of nCount_RNA for each sample
mean.rna <- list()
# create output directories
dir.create(output.metrics, showWarnings = F, recursive = T)
dir.create(output.matrices, showWarnings = F, recursive = T)

# load metadata of scRNAseq samples that will be analyzed
metadata.sample <- read.delim(file = file.path(github.directory, "scRNAseq", "metadata.txt"))
# use names of each samples
all.sorted <- metadata.sample$Sample.name

# create individual seurat objects and subset bad quality cells quantify cell cycle before sample pooling
for (my.sample in all.sorted) {
  all.meta <- metadata.sample[metadata.sample$Sample.name == my.sample, ]
  print(my.sample)
  # read 10x cell ranger filtered matrix
  all.sorted.raw.data[[my.sample]] <- Read10X(file.path(input.directory, my.sample))
  
  # create seurat object filtering genes barely expressed (less than 3 cells) and droplet with minimal amount of UMI (less than 200 features)
  seurat.all.sorted[[my.sample]] <- CreateSeuratObject(counts = all.sorted.raw.data[[my.sample]], project = my.sample, min.cells = 3, 
                                                       min.features = 200)
  
  # Quantify percentage mitochondrial features per cell
  seurat.all.sorted[[my.sample]][["percent.mt"]] <- PercentageFeatureSet(seurat.all.sorted[[my.sample]], pattern = "^mt-")
  # Save QC containing nFeature, nCount and percent.mt displayed in violin Plot
  ggsave(filename = paste0(output.metrics, my.sample, ".RNA.QC.pdf"),
         VlnPlot(seurat.all.sorted[[my.sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3), device = "pdf", width = 7, height = 7)
  
  # quantify mean value of the nCount_RNA to define sample specific threshold for filtering low quality cells/doublets
  mean.rna[[my.sample]] <- mean(seurat.all.sorted[[my.sample]]$nCount_RNA)
  seurat.all.sorted.subset[[my.sample]] <- subset(seurat.all.sorted[[my.sample]], 
                                                  subset = nCount_RNA > 0.4 * mean.rna[[my.sample]] & 
                                                    nCount_RNA < 2.5 * mean.rna[[my.sample]] &
                                                    percent.mt < 8 & percent.mt > 0.05)
  # save metrics post_filtering
  ggsave(filename = paste0(output.metrics, my.sample, ".RNA.QC.post.subset.pdf"),
         VlnPlot(seurat.all.sorted.subset[[my.sample]],
                 features = c("nFeature_RNA",
                              "nCount_RNA",
                              "percent.mt"),
                 ncol = 3), device = "pdf", width = 7, height = 7)
  # NormalizeData, find variableFeatures, and score cell cycle per sample
  seurat.all.sorted.subset[[my.sample]] <- NormalizeData(seurat.all.sorted.subset[[my.sample]])
  seurat.all.sorted.subset[[my.sample]] <- FindVariableFeatures(seurat.all.sorted.subset[[my.sample]],
                                                                nfeatures = 3000, 
                                                                selection.method = "vst")
  seurat.all.sorted.subset[[my.sample]] <- CellCycleScoring( seurat.all.sorted.subset[[my.sample]],
                                                             s.features = s.genes,
                                                             g2m.features = g2m.genes, 
                                                             set.ident = TRUE) 
  # add metadata from the metadata.sample for each samples
  for (current.metadata in colnames(all.meta)) {
    seurat.all.sorted.subset[[my.sample]][[current.metadata]] <- all.meta[[current.metadata]]
  }
}

# combine individual seurat object into a combined one, cell id to use is 
combined.seurat <- merge(seurat.all.sorted.subset[[all.sorted[1]]],
                         y = subsetByNamesOrIndices(seurat.all.sorted.subset, 
                                                    all.sorted[2:length(all.sorted)]),
                         add.cell.ids = all.sorted, 
                         project = "aggregate") 

# If you need to decrease your memory usage:
# rm(list = c("seurat.all.sorted.subset", "all.sorted.raw.data", "seurat.all.sorted"))
# gc()

# set order for ploting seurat object in the same order as the metadata 
for (current.metadata in colnames(metadata.sample)) {
  combined.seurat[[current.metadata]] <- factor(combined.seurat[[]][, current.metadata], 
                                                levels = unique(metadata.sample[match(all.sorted, metadata.sample$Sample.name), current.metadata]))
}

# Visualize number of RNA
VlnPlot(combined.seurat, features = c("nFeature_RNA",
                                      "nCount_RNA",
                                      "percent.mt"),
        group.by = "orig.ident", split.by = "orig.ident", pt.size = 0)
ggsave(filename = paste0(output.metrics, "combined.RNA.QC.post.pdf"), device = "pdf", width = 7, height = 7)

# Normalize Data, find vairable features, scaleData and regress cell cycle score
combined.seurat <- NormalizeData(combined.seurat,
                                 verbose = T)
combined.seurat <- FindVariableFeatures(combined.seurat,
                                        verbose = T,nfeatures = 3000)
combined.seurat <- ScaleData(combined.seurat,
                             verbose = T,
                             vars.to.regress = c("percent.mt", "S.Score",
                                                 "G2M.Score"))
# compute the average gene expression
# and restrict PCA analysis to genes between 20 and 80 percentile over the different samples
average.expression <- rowMeans(AverageExpression(combined.seurat, group.by = "orig.ident")[[1]])
genes.to.exclude <- names(average.expression[average.expression < quantile(average.expression, 0.2) |
                                               average.expression > quantile(average.expression, 0.8)])
used.var.gene <- setdiff(VariableFeatures(combined.seurat),
                         c(genes.to.exclude))

# run PCA analysis
combined.seurat <- RunPCA(combined.seurat, npcs = 50,
                          verbose = FALSE, features = used.var.gene)
# identify number of component to use by analyzing the Elbowplot
ElbowPlot(combined.seurat, ndims = 50)

# Generate UMAP projection for dimensional reduction
combined.seurat <- RunUMAP(combined.seurat, reduction = "pca",
                           n.components = 2L,  
                           dims = 1:15, seed.use = 20)
# build neighbor graph (KNN) using same number of dimension as in UMAP
combined.seurat <- FindNeighbors(combined.seurat, reduction = "pca",
                                 dims = 1:15)
# cluster cells
combined.seurat <- FindClusters(combined.seurat, resolution = 0.6)
# visualize data by time and cell identity
DimPlot(combined.seurat,
        label = T,
        pt.size = 0.1,
        label.size = 7, shuffle = T) +
  DimPlot(combined.seurat,
          label = T,
          pt.size = 0.1,
          label.size = 7,
          group.by = "Time", shuffle = T)

# rename clusters to biologicate fate names
new.cluster.name <- c("Somite A", "Spinal cord A", "NMPs C", "NMPs A",
                      "Early PSM", "Somite B", "Determination front",
                      "Late PSM",
                      "Spinal cord B", "NMPs B", "Endoderm",
                      "Pluripotent/PGC")
names(new.cluster.name) <- levels(combined.seurat)
combined.seurat <- RenameIdents(combined.seurat, new.cluster.name)
combined.seurat[["Fate"]] <- Idents(combined.seurat)
combined.seurat$Fate <- factor(combined.seurat$Fate,levels = c("NMPs A", "NMPs B",
                                                               "NMPs C", "Spinal cord A",
                                                               "Spinal cord B", "Early PSM",
                                                               "Late PSM", "Determination front",
                                                               "Somite A", "Somite B",
                                                               "Endoderm", "Pluripotent/PGC"))

# define color choice and assign to cluster
cluster.color <- c("#DAE3F3", "#8FAADC", "#4472C4",
                   "#2F5597", "#002060",
                   "#FBE5D6", "#F4B183", "#ED7D31",
                   "#C55A11", "#843C0C",
                   "#70AD47", "#7030A0")
names(cluster.color) <- c("NMPs A", "NMPs B", "NMPs C",
                          "Spinal cord A", "Spinal cord B",
                          "Early PSM", "Late PSM", "Determination front",
                          "Somite A", "Somite B",
                          "Endoderm", "Pluripotent/PGC")
combined.seurat$Fate <- factor(combined.seurat$Fate, levels = names(cluster.color))

combined.seurat <- RenameIdents(combined.seurat, cluster.color)
combined.seurat[["Fate.Color"]] <- Idents(combined.seurat)

# For aesthetic reasons we prefer to have Endoderm on top left
endo.cells <- which(combined.seurat[["Fate"]] == "Endoderm")
average.umap1.endo <- mean(Embeddings(object = combined.seurat, reduction = "umap")[endo.cells, 1])
if (average.umap1.endo < 0) {
  combined.seurat[["umap"]]@cell.embeddings[, 1] <- -Embeddings(object = combined.seurat, reduction = "umap")[, 1]
}
average.umap2.endo <- mean(Embeddings(object = combined.seurat, reduction = "umap")[endo.cells, 2])
if (average.umap2.endo < 0) {
  combined.seurat[["umap"]]@cell.embeddings[, 2] <- -Embeddings(object = combined.seurat, reduction = "umap")[, 2]
}

dir.create(file.path(github.directory, "scRNAseq", "plots"), showWarnings = F, recursive = T)
# visualize and save plots of cluster identity and time of gastruloid development
# These plots will be used for Supplementary Figure 6
DimPlot(combined.seurat, group.by = "Fate",
        cols = levels(combined.seurat$Fate.Color),
        pt.size = 0.5, shuffle = T)
ggsave(filename = file.path(github.directory, "scRNAseq", "plots", "Cluster.pdf"), width = 12, height = 8)

DimPlot(combined.seurat,
        label = F,
        pt.size = 0.5, 
        label.size = 4,
        group.by = "Time",
        shuffle = T, cols = grey(level = c(0.8,0.2)))
ggsave(filename = file.path(github.directory, "scRNAseq", "plots", "Time.pdf"), width = 12, height = 8)

# visualize and save plot of desired gene expression
FeaturePlot(combined.seurat, features = c("Hoxd4", "Hoxd8", "Cdx1",
                                          "Hoxd9", "Hoxd10", "Cdx2"),
            ncol = 3, order = T, pt.size = 0.5)
ggsave(filename = file.path(github.directory, "scRNAseq", "plots", "SupFig3C.pdf"), width = 36, height = 16)

# Evaluate the proportion of cells with Cdx1 or Cdx2:
prop.table(table(GetAssayData(object = combined.seurat, assay = "RNA", slot = "counts")["Cdx2",] == 0 &
                   GetAssayData(object = combined.seurat, assay = "RNA", slot = "counts")["Cdx1",] == 0,
                 combined.seurat$Time), 2)
# TRUE = no Cdx1 and no Cdx2
# FALSE = Cdx1 or Cdx2
#              96h       120h
# FALSE 0.97706603 0.56302564
# TRUE  0.02293397 0.43697436

# Save the Seurat object for GEO:
saveRDS(object = combined.seurat, file = file.path(output.matrices, "combined.96h.120h.rds"))
# Save the new metadata:
write.csv(combined.seurat[[]], file = file.path(output.matrices, "metadata.csv"))
# Save the raw counts:
write.csv(as.matrix(GetAssayData(object = combined.seurat, assay = "RNA", slot = "counts")), file = file.path(output.matrices, "raw_counts.csv"))
# Warning: when importing the matrix in R, the colnames ending with '-1' will be substituted by '.1'
