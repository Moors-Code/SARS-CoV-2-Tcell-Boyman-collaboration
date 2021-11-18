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
# Part 3 - Integration of the datasets to generate one dataset that contains all cells
###########################################################################################################################################
# I am following the instructions for integration from:
# https://satijalab.org/seurat/articles/integration_introduction.html
# Create lists containing the datasets
Dataset_list <- list(Set1_sorted,Set2_sorted,Set3_sorted,Set4_sorted,Tcell)


# Select features that are repeatedly variable across datasets for integration
Integreation_features <- SelectIntegrationFeatures(object.list = Dataset_list)


# Perform integration --> Takes long!
Integration_anchors <- FindIntegrationAnchors(object.list = Dataset_list, anchor.features = Integreation_features)


# This command creates an 'integrated' data assay
Integrated <- IntegrateData(anchorset = Integration_anchors)

remove(Dataset_list)
remove(Integration_anchors)
remove(Integreation_features)
remove(Set1_sorted,Set2_sorted,Set3_sorted,Set4_sorted)
remove(Tcell)

Integrated@meta.data$Patient[is.na(Integrated@meta.data$Patient)] <- "Healthy"