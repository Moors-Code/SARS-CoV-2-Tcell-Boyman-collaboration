# Final analysis script of the COVID Tcell scRNAseq project.
# This script produces the scRNAseq plots of the paper.


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
# Part 4 - Dimensional reduction and clustering of the Integrated dataset
###########################################################################################################################################
Integrated <- ScaleData(Integrated)
Integrated <- RunPCA(Integrated, features = VariableFeatures(object = Integrated))
Integrated <- FindNeighbors(Integrated, dims = 1:15,reduction = "pca")
Integrated <- FindClusters(Integrated, resolution = 0.5)
Integrated <- RunUMAP(Integrated, dims = 1:15)

# We enable here the identification of heterogeneously expressed genes across clusters. Commented out for running speed reasons.
#markers <- FindAllMarkers(Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay = "RNA")
#markers.table <- markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
#markers.table <- markers.table[,c(2,5,6,7)]
#write.xlsx(markers.table, 'All_cluster_markers.xlsx')