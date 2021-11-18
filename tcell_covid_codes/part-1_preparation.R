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
# Part 1 - Loading and preparation of the main data
###########################################################################################################################################
# Setting the working directory
setwd("~/NAS/Jan/Experiments and Data/210426 - Exp 033 - CD8 Tcell scRNAseq (Boyman)/cellranger_multi_outputs")

# Loading the data into R
Set1_sorted <- "./multi_Set_1/outs/count/filtered_feature_bc_matrix"
Set2_sorted <- "./multi_Set_2_Sorted/outs/count/filtered_feature_bc_matrix"
Set2_hashed <- "./multi_Set_2_Hashed/outs/count/filtered_feature_bc_matrix"
Set3_sorted <- "./multi_Set_3/outs/count/filtered_feature_bc_matrix"
Set4_hashed <- "./multi_Set_4_Hashed/outs/count/filtered_feature_bc_matrix"
Set4_sorted <- "./multi_Set_4_Sorted/outs/count/filtered_feature_bc_matrix"

Set1_sorted.data <- Read10X(data.dir = Set1_sorted)
Set2_sorted.data <- Read10X(data.dir = Set2_sorted)
Set2_hashed.data <- Read10X(data.dir = Set2_hashed)
Set3_sorted.data <- Read10X(data.dir = Set3_sorted)
Set4_hashed.data <- Read10X(data.dir = Set4_hashed)
Set4_sorted.data <- Read10X(data.dir = Set4_sorted)


# Removing not needed file links
remove(Set1_sorted, Set2_sorted, Set2_hashed, Set3_sorted, Set4_hashed, Set4_sorted)

# Translating HTO IDs into Sample IDs
rownames(Set2_hashed.data$`Antibody Capture`)[10:19] <- c(
  "Patient 1", "Patient 5", "Patient 21", "Patient 25", "Patient 27", "Patient 33", "Patient 35", "Patient 38",
  "Patient 41", "Patient 45"
)
rownames(Set4_hashed.data$`Antibody Capture`)[10:19] <- c(
  "Patient 48", "Patient 50", "Patient 57", "Patient 63", "Patient 68", "Patient 103", "Patient 106", "Patient 112",
  "Patient 120", "Patient 133"
)

# Seurat Objects are created using the Gene expression data.
Set1_sorted <- CreateSeuratObject(counts = Set1_sorted.data$`Gene Expression`)
Set2_sorted <- CreateSeuratObject(counts = Set2_sorted.data$`Gene Expression`)
Set2_hashed <- CreateSeuratObject(counts = Set2_hashed.data$`Gene Expression`)
Set3_sorted <- CreateSeuratObject(counts = Set3_sorted.data$`Gene Expression`)
Set4_hashed <- CreateSeuratObject(counts = Set4_hashed.data$`Gene Expression`)
Set4_sorted <- CreateSeuratObject(counts = Set4_sorted.data$`Gene Expression`)


# One assay is created for the Hashing, one assay is created for the surface protein quantification, one assay is created for the dextramers
# HTO assays (Only for the hashed datasets):
Set2_hashed_HTO_assay <- CreateAssayObject(counts = Set2_hashed.data$`Antibody Capture`[c(10:19), ])
Set4_hashed_HTO_assay <- CreateAssayObject(counts = Set4_hashed.data$`Antibody Capture`[c(11:19), ])

# Protein capture assays (Only for the sorted datasets):
Set1_sorted_Protein_assay <- CreateAssayObject(counts = Set1_sorted.data$`Antibody Capture`[c(8:9), ])
Set2_sorted_Protein_assay <- CreateAssayObject(counts = Set2_sorted.data$`Antibody Capture`[c(8:9), ])
Set3_sorted_Protein_assay <- CreateAssayObject(counts = Set3_sorted.data$`Antibody Capture`[c(8:9), ])
Set4_sorted_Protein_assay <- CreateAssayObject(counts = Set4_sorted.data$`Antibody Capture`[c(8:9), ])

# Dextramer capture assays (Only for the sorted datasets):
Set1_sorted_Dextramer_assay <- CreateAssayObject(counts = Set1_sorted.data$`Antibody Capture`[c(1:7), ])
Set2_sorted_Dextramer_assay <- CreateAssayObject(counts = Set2_sorted.data$`Antibody Capture`[c(1:7), ])
Set3_sorted_Dextramer_assay <- CreateAssayObject(counts = Set3_sorted.data$`Antibody Capture`[c(1:7), ])
Set4_sorted_Dextramer_assay <- CreateAssayObject(counts = Set4_sorted.data$`Antibody Capture`[c(1:7), ])

# Now the assays are added to the previously created Seurat objects
Set2_hashed[["HTO"]] <- Set2_hashed_HTO_assay
Set4_hashed[["HTO"]] <- Set4_hashed_HTO_assay

Set1_sorted[["Protein"]] <- Set1_sorted_Protein_assay
Set2_sorted[["Protein"]] <- Set2_sorted_Protein_assay
Set3_sorted[["Protein"]] <- Set3_sorted_Protein_assay
Set4_sorted[["Protein"]] <- Set4_sorted_Protein_assay

Set1_sorted[["Dextramer"]] <- Set1_sorted_Dextramer_assay
Set2_sorted[["Dextramer"]] <- Set2_sorted_Dextramer_assay
Set3_sorted[["Dextramer"]] <- Set3_sorted_Dextramer_assay
Set4_sorted[["Dextramer"]] <- Set4_sorted_Dextramer_assay

remove(
  Set1_sorted_Dextramer_assay, Set1_sorted_Protein_assay, Set2_hashed_HTO_assay, Set2_sorted_Dextramer_assay,
  Set2_sorted_Protein_assay, Set3_sorted_Dextramer_assay, Set3_sorted_Protein_assay, Set4_hashed_HTO_assay, Set4_sorted_Dextramer_assay, Set4_sorted_Protein_assay
)

remove(Set1_sorted.data, Set2_sorted.data, Set2_hashed.data, Set3_sorted.data, Set4_hashed.data, Set4_sorted.data)

# QC and selecting cells for further analysis
Set1_sorted[["percent.mt"]] <- PercentageFeatureSet(Set1_sorted, pattern = "^MT-")
Set2_sorted[["percent.mt"]] <- PercentageFeatureSet(Set2_sorted, pattern = "^MT-")
Set2_hashed[["percent.mt"]] <- PercentageFeatureSet(Set2_hashed, pattern = "^MT-")
Set3_sorted[["percent.mt"]] <- PercentageFeatureSet(Set3_sorted, pattern = "^MT-")
Set4_hashed[["percent.mt"]] <- PercentageFeatureSet(Set4_hashed, pattern = "^MT-")
Set4_sorted[["percent.mt"]] <- PercentageFeatureSet(Set4_sorted, pattern = "^MT-")

Set1_sorted <- subset(Set1_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set2_sorted <- subset(Set2_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set2_hashed <- subset(Set2_hashed, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set3_sorted <- subset(Set3_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set4_hashed <- subset(Set4_hashed, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Set4_sorted <- subset(Set4_sorted, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# How many cells are left in each of the datasets?
# nrow(Set1_sorted@meta.data)
# nrow(Set2_sorted@meta.data)
# nrow(Set2_hashed@meta.data)
# nrow(Set3_sorted@meta.data)
# nrow(Set4_hashed@meta.data)
# nrow(Set4_sorted@meta.data)

# Now, normalization of the Gene expression data and of the assay data is done.
Set1_sorted <- NormalizeData(Set1_sorted)
Set2_sorted <- NormalizeData(Set2_sorted)
Set2_hashed <- NormalizeData(Set2_hashed)
Set3_sorted <- NormalizeData(Set3_sorted)
Set4_hashed <- NormalizeData(Set4_hashed)
Set4_sorted <- NormalizeData(Set4_sorted)

Set2_hashed <- NormalizeData(Set2_hashed, assay = "HTO", normalization.method = "CLR")
Set4_hashed <- NormalizeData(Set4_hashed, assay = "HTO", normalization.method = "CLR")

Set1_sorted <- NormalizeData(Set1_sorted, assay = "Protein", normalization.method = "CLR", margin = 2)
Set2_sorted <- NormalizeData(Set2_sorted, assay = "Protein", normalization.method = "CLR", margin = 2)
Set3_sorted <- NormalizeData(Set3_sorted, assay = "Protein", normalization.method = "CLR", margin = 2)
Set4_sorted <- NormalizeData(Set4_sorted, assay = "Protein", normalization.method = "CLR", margin = 2)

# Demultiplexing the HTO data
Set2_hashed <- HTODemux(Set2_hashed, assay = "HTO", positive.quantile = 0.99)
Set4_hashed <- HTODemux(Set4_hashed, assay = "HTO", positive.quantile = 0.99)

# table(Set2_hashed@meta.data$hash.ID)
# table(Set4_hashed@meta.data$hash.ID)

# df <- as.data.frame(table(Set2_hashed@meta.data$hash.ID))
# df2 <- as.data.frame(table(Set4_hashed@meta.data$hash.ID))
# ggplot(df, aes(Var1, Freq)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   ggtitle("Set 2 hashed") +
#   theme(plot.title = element_text(hjust = 0.5))
# ggplot(df2, aes(Var1, Freq)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   ggtitle("Set 4 hashed") +
#   theme(plot.title = element_text(hjust = 0.5))
# remove(df, df2)

# Adding patient identities to cells using the souporcell output
# 1) Load data into
SNPs.1_2 <- read.csv("snp_demux_set_1_2_merged/clusters.tsv", sep = "\t")
SNPs.1_2 <- SNPs.1_2[, c(1:3)]
# table(SNPs.1_2$status)
SNPs.1_2 <- SNPs.1_2[SNPs.1_2$status == "singlet", ]

SNPs.3_4 <- read.csv("snp_demux_set_3_4_merged/clusters.tsv", sep = "\t")
SNPs.3_4 <- SNPs.3_4[, c(1:3)]
# table(SNPs.3_4$status)
SNPs.3_4 <- SNPs.3_4[SNPs.3_4$status == "singlet", ]

# 2) merge hashed singlets with souporcell output
hashed.df2 <- Set2_hashed@meta.data
hashed.df4 <- Set4_hashed@meta.data

rownames(hashed.df2) <- paste0(rownames(hashed.df2), "_Set_2_Hashed")
rownames(hashed.df4) <- paste0(rownames(hashed.df4), "_Set_4_Hashed")

hashed.df2 <- hashed.df2[hashed.df2$HTO_classification.global == "Singlet", ]
hashed.df4 <- hashed.df4[hashed.df4$HTO_classification.global == "Singlet", ]

hashed.df2$barcode <- rownames(hashed.df2)
hashed.df4$barcode <- rownames(hashed.df4)

SNPs.1.2.merged <- merge(SNPs.1_2, hashed.df2, by = "barcode", all.x = TRUE)
SNPs.3.4.merged <- merge(SNPs.3_4, hashed.df4, by = "barcode", all.x = TRUE)


SNPs.1.2.merged <- na.omit(SNPs.1.2.merged)
SNPs.3.4.merged <- na.omit(SNPs.3.4.merged)

# 3) Test if same singlets are called to be same souporcell clusters
# coul <- brewer.pal(10, "Set3")
# 
# tbl <- with(SNPs.1.2.merged, table(assignment, hash.ID))
# barplot(tbl, beside = TRUE, legend = TRUE, col = coul)
# 
# tbl4 <- with(SNPs.3.4.merged, table(assignment, hash.ID))
# barplot(tbl4, beside = TRUE, legend = TRUE, col = coul)
# 
# remove(coul, tbl, tbl4)

# 4) Add patient info to datasets, for those patients where imputation can be done
SNPs.1_2$Patient <- "unknown"
SNPs.3_4$Patient <- "unknown"

SNPs.1_2$Patient[SNPs.1_2$assignment == 0] <- "Patient 1"
SNPs.1_2$Patient[SNPs.1_2$assignment == 1] <- "Patient 33"
SNPs.1_2$Patient[SNPs.1_2$assignment == 2] <- "Patient 21"
SNPs.1_2$Patient[SNPs.1_2$assignment == 5] <- "Patient 25"
SNPs.1_2$Patient[SNPs.1_2$assignment == 6] <- "Patient 38"
SNPs.1_2$Patient[SNPs.1_2$assignment == 7] <- "Patient 27"
SNPs.1_2$Patient[SNPs.1_2$assignment == 9] <- "Patient 45"

SNPs.3_4$Patient[SNPs.3_4$assignment == 0] <- "Patient 103"
SNPs.3_4$Patient[SNPs.3_4$assignment == 1] <- "Patient 63"
SNPs.3_4$Patient[SNPs.3_4$assignment == 2] <- "Patient 112"
SNPs.3_4$Patient[SNPs.3_4$assignment == 3] <- "Patient 120"
SNPs.3_4$Patient[SNPs.3_4$assignment == 4] <- "Patient 50"
SNPs.3_4$Patient[SNPs.3_4$assignment == 5] <- "Patient 57"
SNPs.3_4$Patient[SNPs.3_4$assignment == 8] <- "Patient 106"
SNPs.3_4$Patient[SNPs.3_4$assignment == 9] <- "Patient 48"

# Making vectors containing rownames (part of preparation for the merging by rowname, which follows after)
rownames.Set1_sorted <- rownames(Set1_sorted@meta.data)
rownames.Set2_sorted <- rownames(Set2_sorted@meta.data)
rownames.Set3_sorted <- rownames(Set3_sorted@meta.data)
rownames.Set4_sorted <- rownames(Set4_sorted@meta.data)

rownames(Set1_sorted@meta.data) <- paste0(rownames(Set1_sorted@meta.data), "_Set_1")
rownames(Set2_sorted@meta.data) <- paste0(rownames(Set2_sorted@meta.data), "_Set_2_Sorted")
rownames(Set3_sorted@meta.data) <- paste0(rownames(Set3_sorted@meta.data), "_Set_3")
rownames(Set4_sorted@meta.data) <- paste0(rownames(Set4_sorted@meta.data), "_Set_4_Sorted")

rownames(SNPs.1_2) <- SNPs.1_2$barcode
rownames(SNPs.3_4) <- SNPs.3_4$barcode

# Do the merge and replace NA values with 'unknown'
Set1_sorted@meta.data <- merge(Set1_sorted@meta.data, SNPs.1_2, by = 0, all.x = TRUE)
Set2_sorted@meta.data <- merge(Set2_sorted@meta.data, SNPs.1_2, by = 0, all.x = TRUE)
Set3_sorted@meta.data <- merge(Set3_sorted@meta.data, SNPs.3_4, by = 0, all.x = TRUE)
Set4_sorted@meta.data <- merge(Set4_sorted@meta.data, SNPs.3_4, by = 0, all.x = TRUE)

# nrow(Set1_sorted@meta.data)
# table(Set1_sorted@meta.data$status, useNA = "always")
# table(Set2_sorted@meta.data$status, useNA = "always")
# table(Set3_sorted@meta.data$status, useNA = "always")
# table(Set4_sorted@meta.data$status, useNA = "always")

Set1_sorted@meta.data$Patient[is.na(Set1_sorted@meta.data$Patient)] <- "unknown"
Set2_sorted@meta.data$Patient[is.na(Set2_sorted@meta.data$Patient)] <- "unknown"
Set3_sorted@meta.data$Patient[is.na(Set3_sorted@meta.data$Patient)] <- "unknown"
Set4_sorted@meta.data$Patient[is.na(Set4_sorted@meta.data$Patient)] <- "unknown"

Set1_sorted@meta.data$status[is.na(Set1_sorted@meta.data$status)] <- "unknown"
Set2_sorted@meta.data$status[is.na(Set2_sorted@meta.data$status)] <- "unknown"
Set3_sorted@meta.data$status[is.na(Set3_sorted@meta.data$status)] <- "unknown"
Set4_sorted@meta.data$status[is.na(Set4_sorted@meta.data$status)] <- "unknown"

Set1_sorted@meta.data$assignment[is.na(Set1_sorted@meta.data$assignment)] <- "unknown"
Set2_sorted@meta.data$assignment[is.na(Set2_sorted@meta.data$assignment)] <- "unknown"
Set3_sorted@meta.data$assignment[is.na(Set3_sorted@meta.data$assignment)] <- "unknown"
Set4_sorted@meta.data$assignment[is.na(Set4_sorted@meta.data$assignment)] <- "unknown"

rownames(Set1_sorted@meta.data) <- rownames.Set1_sorted
rownames(Set2_sorted@meta.data) <- rownames.Set2_sorted
rownames(Set3_sorted@meta.data) <- rownames.Set3_sorted
rownames(Set4_sorted@meta.data) <- rownames.Set4_sorted

remove(rownames.Set1_sorted, rownames.Set2_sorted, rownames.Set3_sorted, rownames.Set4_sorted)
remove(hashed.df2, hashed.df4, SNPs.1_2, SNPs.3_4, SNPs.1.2.merged, SNPs.3.4.merged)
remove(Set2_hashed, Set4_hashed)

# Adding suffix to SNP assignment to keep info from which Souporcell run it comes
Set1_sorted@meta.data$assignment <- paste0(Set1_sorted@meta.data$assignment, "_1")
Set2_sorted@meta.data$assignment <- paste0(Set2_sorted@meta.data$assignment, "_1")
Set3_sorted@meta.data$assignment <- paste0(Set3_sorted@meta.data$assignment, "_2")
Set4_sorted@meta.data$assignment <- paste0(Set4_sorted@meta.data$assignment, "_2")


# Perform cell type annotation (matching cell identities previously defined using Azimuth)
Set1_sorted.predictions_1 <- read.delim("./Azimuth results/Set1_sorted/azimuth_pred.tsv", row.names = 1)
Set2_sorted.predictions_1 <- read.delim("./Azimuth results/Set2_sorted/azimuth_pred.tsv", row.names = 1)
Set3_sorted.predictions_1 <- read.delim("./Azimuth results/Set3_sorted/azimuth_pred.tsv", row.names = 1)
Set4_sorted.predictions_1 <- read.delim("./Azimuth results/Set4_sorted/azimuth_pred.tsv", row.names = 1)

Set1_sorted <- AddMetaData(
  object = Set1_sorted,
  metadata = Set1_sorted.predictions_1
)

Set2_sorted <- AddMetaData(
  object = Set2_sorted,
  metadata = Set2_sorted.predictions_1
)

Set3_sorted <- AddMetaData(
  object = Set3_sorted,
  metadata = Set3_sorted.predictions_1
)

Set4_sorted <- AddMetaData(
  object = Set4_sorted,
  metadata = Set4_sorted.predictions_1
)

remove(Set1_sorted.predictions_1, Set2_sorted.predictions_1, Set3_sorted.predictions_1, Set4_sorted.predictions_1)

# Addition of the TCR data
Set1_sorted_TCR <- read.csv("./multi_Set_1/outs/vdj_t/filtered_contig_annotations.csv")
Set2_sorted_TCR <- read.csv("./multi_Set_2_Sorted/outs/vdj_t/filtered_contig_annotations.csv")
Set3_sorted_TCR <- read.csv("./multi_Set_3/outs/vdj_t/filtered_contig_annotations.csv")
Set4_sorted_TCR <- read.csv("./multi_Set_4_Sorted/outs/vdj_t/filtered_contig_annotations.csv")

TCRcombined1 <- combineTCR(Set1_sorted_TCR, samples = "A1", ID = "A1", cells = "T-AB", filterMulti = T)
TCRcombined2 <- combineTCR(Set2_sorted_TCR, samples = "A1", ID = "A1", cells = "T-AB", filterMulti = T)
TCRcombined3 <- combineTCR(Set3_sorted_TCR, samples = "A1", ID = "A1", cells = "T-AB", filterMulti = T)
TCRcombined4 <- combineTCR(Set4_sorted_TCR, samples = "A1", ID = "A1", cells = "T-AB", filterMulti = T)

TCRcombined1$A1_A1$barcode <- substring(TCRcombined1$A1_A1$barcode, 7)
TCRcombined2$A1_A1$barcode <- substring(TCRcombined2$A1_A1$barcode, 7)
TCRcombined3$A1_A1$barcode <- substring(TCRcombined3$A1_A1$barcode, 7)
TCRcombined4$A1_A1$barcode <- substring(TCRcombined4$A1_A1$barcode, 7)

Set1_sorted <- combineExpression(TCRcombined1, Set1_sorted, cloneCall = "aa", cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500))
Set2_sorted <- combineExpression(TCRcombined2, Set2_sorted, cloneCall = "aa", cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500))
Set3_sorted <- combineExpression(TCRcombined3, Set3_sorted, cloneCall = "aa", cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500))
Set4_sorted <- combineExpression(TCRcombined4, Set4_sorted, cloneCall = "aa", cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500))

remove(Set1_sorted_TCR, Set2_sorted_TCR, Set3_sorted_TCR, Set4_sorted_TCR)
remove(TCRcombined1, TCRcombined2, TCRcombined3, TCRcombined4)

# Are there same clones based on "CTaa" in both timepoints? - These will be defined as common clones.
common_Set1.2 <- intersect(Set1_sorted@meta.data$CTaa, Set2_sorted@meta.data$CTaa)
common_Set3.4 <- intersect(Set3_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)
common_Set1.2 <- common_Set1.2[!is.na(common_Set1.2)]
common_Set3.4 <- common_Set3.4[!is.na(common_Set3.4)]

# # Are there overlaps of clones between datasets?
# common_Set1.3 <- intersect(Set1_sorted@meta.data$CTaa, Set3_sorted@meta.data$CTaa)
# common_Set1.4 <- intersect(Set1_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)
# common_Set2.3 <- intersect(Set2_sorted@meta.data$CTaa, Set3_sorted@meta.data$CTaa)
# common_Set2.4 <- intersect(Set2_sorted@meta.data$CTaa, Set4_sorted@meta.data$CTaa)
# 
# common_Set1.3 <- common_Set1.3[!is.na(common_Set1.3)]
# common_Set1.4 <- common_Set1.4[!is.na(common_Set1.4)]
# common_Set2.3 <- common_Set2.3[!is.na(common_Set2.3)]
# common_Set2.4 <- common_Set2.4[!is.na(common_Set2.4)]
# 
# Overlaps <- c(common_Set1.3, common_Set1.4, common_Set2.3)
# Overlaps <- unique(Overlaps)
remove(Overlaps)
remove(common_Set1.3, common_Set1.4, common_Set2.3, common_Set2.4)

# The data.frames that are created here, include from each set the cells which have assigned clonotypes that are found
# in the respective first and second timepoint.
Set1_commons <- filter(Set1_sorted@meta.data, grepl(paste(common_Set1.2, collapse = "|"), CTaa))
Set2_commons <- filter(Set2_sorted@meta.data, grepl(paste(common_Set1.2, collapse = "|"), CTaa))
Set3_commons <- filter(Set3_sorted@meta.data, grepl(paste(common_Set3.4, collapse = "|"), CTaa))
Set4_commons <- filter(Set4_sorted@meta.data, grepl(paste(common_Set3.4, collapse = "|"), CTaa))

remove(common_Set1.2, common_Set3.4)

# Adding a column to the metadata that specifies whether a cell is among the once that have assigned clonotypes that are
# found in the respective first and second timepoint (- basically adding a "yes" to all cells that appear in
# Set1_commons 1-4).

Set1_sorted@meta.data$common.clone <- "no"
Set2_sorted@meta.data$common.clone <- "no"
Set3_sorted@meta.data$common.clone <- "no"
Set4_sorted@meta.data$common.clone <- "no"

Set1_sorted@meta.data[rownames(Set1_commons), 23] <- "yes"
Set2_sorted@meta.data[rownames(Set2_commons), 23] <- "yes"
Set3_sorted@meta.data[rownames(Set3_commons), 23] <- "yes"
Set4_sorted@meta.data[rownames(Set4_commons), 23] <- "yes"

remove(Set1_commons, Set2_commons, Set3_commons, Set4_commons)


#############################################
#Part 1.1 - Prepare dataset integration
#############################################
# Perform FindVariableFeatures
Set1_sorted <- FindVariableFeatures(Set1_sorted, selection.method = "vst", nfeatures = 2000)
Set2_sorted <- FindVariableFeatures(Set2_sorted, selection.method = "vst", nfeatures = 2000)
Set3_sorted <- FindVariableFeatures(Set3_sorted, selection.method = "vst", nfeatures = 2000)
Set4_sorted <- FindVariableFeatures(Set4_sorted, selection.method = "vst", nfeatures = 2000)

# Rename cells to keep dataset info
Set1_sorted <- RenameCells(object = Set1_sorted, add.cell.id = "Set1")
Set2_sorted <- RenameCells(object = Set2_sorted, add.cell.id = "Set2")
Set3_sorted <- RenameCells(object = Set3_sorted, add.cell.id = "Set3")
Set4_sorted <- RenameCells(object = Set4_sorted, add.cell.id = "Set4")

# Add columns to keep dataset info
Set1_sorted$Dataset <- "Set1"
Set2_sorted$Dataset <- "Set2"
Set3_sorted$Dataset <- "Set3"
Set4_sorted$Dataset <- "Set4"

Set1_sorted$Timepoint <- "Acute"
Set2_sorted$Timepoint <- "Recovered"
Set3_sorted$Timepoint <- "Acute"
Set4_sorted$Timepoint <- "Recovered"
