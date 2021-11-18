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
# Part 11 - DGEA
###########################################################################################################################################
# DGEA Persistent vs non persistent cells
persistent.vs.non.persistent.markers <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents, assay = "RNA"))
persistent.vs.non.persistent.markers <- persistent.vs.non.persistent.markers[,c(2,5)]
persistent.vs.non.persistent.markers <- persistent.vs.non.persistent.markers[1:170,]
write.xlsx(persistent.vs.non.persistent.markers, 'persistent.vs.non.persistent.markers.xlsx',row.names = TRUE)

# DGEA for selected genes
features <- c('IL7R','GZMB','CX3CR1','IL2RG','NR4A1','ZAP70','STAT3')
persistent.vs.non.persistent.markers.selected.genes <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents,
                                                                                 logfc.threshold = 0, assay = "RNA",min.pct = 0,features = features))
persistent.vs.non.persistent.markers.selected.genes <- persistent.vs.non.persistent.markers.selected.genes[,c(2,5)]
write.xlsx(persistent.vs.non.persistent.markers.selected.genes, 'persistent.vs.non.persistent.markers.selected.genes.xlsx',row.names = TRUE)

# DGEA for TotalSeq proteins
features_found.in.Protein <- c('CD45RA','CCR7.1')
persistent.vs.non.persistent.proteins <- as.data.frame(FindMarkers(Integrated_NA_filtered ,ident.1 = Persistents, ident.2 = Non_persistents,
                                                                   logfc.threshold = 0, assay = "Protein",min.pct = 0,features = features_found.in.Protein))
persistent.vs.non.persistent.proteins <- persistent.vs.non.persistent.proteins[,c(2,5)]
write.xlsx(persistent.vs.non.persistent.proteins, 'persistent.vs.non.persistent.proteins.xlsx',row.names = TRUE)

remove(persistent.vs.non.persistent.markers,persistent.vs.non.persistent.markers.selected.genes,persistent.vs.non.persistent.proteins,features,features_found.in.Protein,Non_persistents)
