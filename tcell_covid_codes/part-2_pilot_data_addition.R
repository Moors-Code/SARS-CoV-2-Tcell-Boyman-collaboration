# Loading packages
library(openxlsx)
library(presto)
library(msigdbr)
library(fgsea)
library(Seurat)
library(scRepertoire)
library(umap)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(dplyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(VennDiagram)
library(patchwork)
library(ggpubr)
library(pals)
library(ggvenn)
library(eulerr)
options(ggrepel.max.overlaps = Inf)

# Clearing the environment
rm(list = ls())
options(scipen = 999)
options(scipen = 0)


###########################################################################################################################################
# Part 2 - Loading and preparation of the pilot data (take healthy patient cells only)
###########################################################################################################################################
# Creating links to the data.
data_dir <- "./Pilot data/filtered_feature_bc_matrix"
SNPs_clusters <- read.csv("./Pilot data/souporcell_5samples.TCR_Boyman_mar2021/snp_demux_Tcells_Boyman_allCustom/clusters.tsv", sep = "\t")

# Loading the data into R
Tcells <- Read10X(data.dir = data_dir)
remove(data_dir)

# Here the Seurat object is created and the counts are normalized
Tcell <- CreateSeuratObject(counts = Tcells$`Gene Expression`)
Tcell <- NormalizeData(Tcell)

# Here, the HTO data and the TotalSeq data and the dCODE data are added as independent assays.
Tcell[["HTO"]] <- CreateAssayObject(counts = Tcells$Custom[c(1:10), ])
Tcell[["dCODEs"]] <- CreateAssayObject(counts = Tcells$Custom[c(12:14), ])
Tcell[["SurfaceProtein"]] <- CreateAssayObject(counts = as.matrix(Tcells$Custom)[11, , drop = FALSE])

# Also, the these assay objects are normalized.
# Here, it is not clear to me yet, for which assays this is necessary.
Tcell <- NormalizeData(Tcell, assay = "HTO", normalization.method = "CLR")
Tcell <- NormalizeData(Tcell, assay = "SurfaceProtein", normalization.method = "CLR")

# Data filtering
Tcell[["percent.mt"]] <- PercentageFeatureSet(Tcell, pattern = "^MT-")

# This is the step where the filtering by mitochondrial genes and nFeature counts happens.
# For this dataset, I have decided to go with a filtering step that exludes cells with >10% mitochondrial genes.
Tcell <- subset(Tcell, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# In this step, the cells are demultiplexed based on HTO enrichment.
Tcell <- HTODemux(Tcell, assay = "HTO", positive.quantile = 0.99)
table(Tcell@meta.data$hash.ID)

# Now the vdj data is added.
# First, the clonotypedata output file from cellranger is loaded into R
vdj <- read.csv("./Pilot data/filtered_contig_annotations.csv")

# As the output of CellRanger are quantifications of both chains,
# the next step is to create a single list object with the by cell barcode.
vdjcombined1 <- combineTCR(vdj, samples = "A1", ID = "A1")

# Then, the clonotypedata is paired with the seurat object that we already have.
# Importantly, the major requirement for the attachment is matching contig cell
# barcodes and barcodes in the row names of the meta data of the seurat object.
# If these do not match, the attachment will fail.
vdjcombined1$A1_A1$barcode <- substring(vdjcombined1$A1_A1$barcode, 7)
Tcell <- combineExpression(vdjcombined1, Tcell, cloneCall = "CTaa")

# Now I combine the HTO demultiplexed data with the SNP demultiplexed data.
shortSNPs_clusters <- SNPs_clusters[, c(1:3)]
Tcell@meta.data$barcode <- rownames(Tcell@meta.data)
Tcell@meta.data <- merge(Tcell@meta.data, shortSNPs_clusters, by = "barcode")
rownames(Tcell@meta.data) <- Tcell@meta.data$barcode
Tcell@meta.data

df <- as.data.frame(count(Tcell@meta.data, HTO_classification.global))
plot1 <- ggplot(df, aes(HTO_classification.global, n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(label = n), size = 3, hjust = 0.5, vjust = -1, position = "stack")
plot1
df2 <- as.data.frame(count(Tcell@meta.data, status))
plot2 <- ggplot(df2, aes(status, n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(label = n), size = 3, hjust = 0.5, vjust = -1, position = "stack")
plot2
plot1 + plot2

Tcell2 <- subset(Tcell, subset = HTO_classification.global == "Singlet")
Tcell2 <- subset(Tcell2, subset = status == "singlet")
nrow(Tcell2@meta.data)
head(Tcell2@meta.data)
tbl <- with(Tcell2@meta.data, table(assignment, HTO_maxID))
tbl
barplot(tbl, beside = TRUE, legend = TRUE)

# The plot above shows that the SNPs and the HTO data resemble each other very well. I am going to keep all healthy donor cells,
# which are SNPs clusters 1,2,3,4
keep.clusters <- c(1, 2, 3, 4)
Tcell <- subset(Tcell, subset = assignment %in% keep.clusters)
Tcell@meta.data$assignment <- "Healthy"

# FindVariableFeatures is performed.
Tcell <- FindVariableFeatures(Tcell, selection.method = "mean.var.plot")

# Preparing for integration
Tcell <- RenameCells(object = Tcell, add.cell.id = "Healthy")
Tcell$Dataset <- "Healthy"
Tcell$Timepoint <- "Healthy"

# Removal of unnessecary files
remove(df, df2, plot1, plot2, shortSNPs_clusters, SNPs_clusters, Tcell2, Tcells, vdj, vdjcombined1, keep.clusters, tbl)

