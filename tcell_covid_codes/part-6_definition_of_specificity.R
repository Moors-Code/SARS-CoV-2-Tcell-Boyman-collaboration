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
# Part 6 - Dextramer part --> Define whether a cell is specific for a Dextramer
###########################################################################################################################################
# The first task is to remove the cells for which we don't have a TCR now.
Integrated_NA_filtered <- Integrated
Integrated_NA_filtered@meta.data$CTaa.known <- "yes"
Integrated_NA_filtered@meta.data$CTaa.known[is.na(Integrated_NA_filtered@meta.data$CTaa)] <- "no"
Integrated_NA_filtered <- subset(Integrated_NA_filtered, subset= CTaa.known == "yes")


Dextramer_counts <- as.data.frame(Integrated_NA_filtered@assays$Dextramer@counts)
Dextramer_counts_transposed <- data.table::transpose(Dextramer_counts)
rownames(Dextramer_counts_transposed) <- colnames(Dextramer_counts)
colnames(Dextramer_counts_transposed) <- rownames(Dextramer_counts)
Dextramer_counts <- Dextramer_counts_transposed
remove(Dextramer_counts_transposed)

# UMI count needs to be higher than 10 AND more than 5x higher than its negative control.
# Calculating ratios for each Dextramer vs. its negative control
Dextramer_counts$FTSDYYQLY.corrected <- Dextramer_counts$FTSDYYQLY/Dextramer_counts$`negative-Ctrl-A01-01`
Dextramer_counts$TTDPSFLGRY.corrected <- Dextramer_counts$TTDPSFLGRY/Dextramer_counts$`negative-Ctrl-A01-01`
Dextramer_counts$ATEGALNTPK.corrected <- Dextramer_counts$ATEGALNTPK/Dextramer_counts$`negative-Ctrl-Universal`
Dextramer_counts$KTFPPTEPK.corrected <- Dextramer_counts$KTFPPTEPK/Dextramer_counts$`negative-Ctrl-Universal`

# Classifying if a cell is Positive for a specific Dextramer
Dextramer_counts$FTSDYYQLY.classification <- "Negative"
Dextramer_counts$TTDPSFLGRY.classification <- "Negative"
Dextramer_counts$ATEGALNTPK.classification <- "Negative"
Dextramer_counts$KTFPPTEPK.classification <- "Negative"

Dextramer_counts$FTSDYYQLY.classification[Dextramer_counts$FTSDYYQLY>10 & (Dextramer_counts$FTSDYYQLY.corrected>5 | Dextramer_counts$FTSDYYQLY.corrected==Inf)] <- "Positive"
Dextramer_counts$TTDPSFLGRY.classification[Dextramer_counts$TTDPSFLGRY>10 & (Dextramer_counts$TTDPSFLGRY.corrected>5 | Dextramer_counts$TTDPSFLGRY.corrected==Inf)] <- "Positive"
Dextramer_counts$ATEGALNTPK.classification[Dextramer_counts$ATEGALNTPK>10 & (Dextramer_counts$ATEGALNTPK.corrected>5 | Dextramer_counts$ATEGALNTPK.corrected==Inf)] <- "Positive"
Dextramer_counts$KTFPPTEPK.classification[Dextramer_counts$KTFPPTEPK>10 & (Dextramer_counts$KTFPPTEPK.corrected>5 | Dextramer_counts$KTFPPTEPK.corrected==Inf)] <- "Positive"

# For how many Dextramers is a cell Positive?
Dextramer_counts$positive.for.n.dextramers <- rowSums(Dextramer_counts == "Positive") 

# Add the information whether a cell is Dextramer specific to the datasets
# colnames(Integrated_NA_filtered@meta.data)
# colnames(Dextramer_counts)
Integrated_NA_filtered@meta.data[42:46] <- Dextramer_counts[12:16]

# What about double Positive cells? --> Kick them out --> Then the frequency column has to be updated!
Integrated_NA_filtered <- subset(Integrated_NA_filtered, subset = positive.for.n.dextramers<2)
Integrated_NA_filtered@meta.data$Row.names <- rownames(Integrated_NA_filtered@meta.data)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data %>% group_by(Dataset, CTaa) %>% add_count(CTaa,name = "updated.frequency") %>% ungroup()
Integrated_NA_filtered@meta.data$Frequency <- Integrated_NA_filtered@meta.data$updated.frequency
Integrated_NA_filtered@meta.data <- as.data.frame(Integrated_NA_filtered@meta.data)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names

# Scatter plot examples
# FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "negative-Ctrl-A01-01", pt.size = 1,group.by = 'ATEGALNTPK.classification')+
#   ggtitle("ATEGALNTPK classification of cells")+xlab("Dex-ATEGALNTPK UMI counts")+ylab("Dex-control UMI counts")+
#   labs(color = "CoV2-Dex classification")
# FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "negative-Ctrl-A01-01", pt.size = 1,group.by = 'ATEGALNTPK.classification')+xlim(0,30)+
#   ggtitle("ATEGALNTPK classification of cells")+xlab("Dex-ATEGALNTPK UMI counts")+ylab("Dex-control UMI counts")+
#   labs(color = "CoV2-Dex classification")
# FeatureScatter(Integrated_NA_filtered, feature1 = "ATEGALNTPK", feature2 = "KTFPPTEPK", pt.size = 1,group.by = 'positive.for.n.dextramers')+
#   labs(color = "CoV2-Dex classification")+ggtitle("ATEGALNTPK vs. KTFPPTEPK counts of cells")+
#   scale_color_manual(labels = c("CoV2-Dex-", "CoV2-Dex+"), values = c("#F8766D", "#00BFC4"))

remove(Dextramer_counts)

# What about TCRs that appear to be specific for more than one peptide?
# First, how many of them do we have?
# colnames(Integrated_NA_filtered@meta.data)
# Investigation <- Integrated_NA_filtered@meta.data[,c(19,42:46)]
# Investigation <- Investigation[Investigation$FTSDYYQLY.classification=="Positive" | 
#                                  Investigation$TTDPSFLGRY.classification=="Positive" | 
#                                  Investigation$ATEGALNTPK.classification=="Positive"|
#                                  Investigation$KTFPPTEPK.classification=="Positive",]
# Investigation <- unique(Investigation)
# Investigation$CTaa[duplicated(Investigation$CTaa)] 
# remove(Investigation)

# This tells us that TCRs "CAVPKNTGNQFYF_CAISVGNEQFF" and "CALSEGNNDMRF_CATSRGLASTDTQYF" are both Positive for more than one peptide.
# Which peptides are that?
# colnames(Integrated_NA_filtered@meta.data)
# unique(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF",c(42:45)])
# unique(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF",c(42:45)])

# We see that both are Positive for TTDPSFLGRY and for ATEGALNTPK
# I am using Featurescatter plots to investigate further.
# df <- subset(Integrated_NA_filtered, subset = CTaa =="CAVPKNTGNQFYF_CAISVGNEQFF")
# FeatureScatter(df, feature1 = "TTDPSFLGRY", feature2 = "ATEGALNTPK", pt.size = 2,group.by = 'positive.for.n.dextramers')+ggtitle("TCR CAVPKNTGNQFYF_CAISVGNEQFF specificity")
# 
# df <- subset(Integrated_NA_filtered, subset = CTaa =="CALSEGNNDMRF_CATSRGLASTDTQYF")
# FeatureScatter(df, feature1 = "TTDPSFLGRY", feature2 = "ATEGALNTPK", pt.size = 2,group.by = 'positive.for.n.dextramers')+ggtitle("TCR CALSEGNNDMRF_CATSRGLASTDTQYF specificity")
# remove(df)

# From this analyis we decided to manually curate all cells with TCR CAVPKNTGNQFYF_CAISVGNEQFF as being positive for no peptide
# Furthermore, we decided to manually curate cells with TCR CALSEGNNDMRF_CATSRGLASTDTQYF as being positive for peptide ATEGALNTPK
Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- 0
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- "Negative"
Integrated_NA_filtered@meta.data$ATEGALNTPK.classification[Integrated_NA_filtered@meta.data$CTaa=="CAVPKNTGNQFYF_CAISVGNEQFF"] <- "Negative"
Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF" & Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification=="Positive"] <- 0
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$CTaa=="CALSEGNNDMRF_CATSRGLASTDTQYF"] <- "Negative"

# Now the investigation from above should not find anything anymore.
# colnames(Integrated_NA_filtered@meta.data)
# Investigation <- Integrated_NA_filtered@meta.data[,c(19,42:46)]
# Investigation <- Investigation[Investigation$FTSDYYQLY.classification=="Positive"|Investigation$TTDPSFLGRY.classification=="Positive"|Investigation$ATEGALNTPK.classification=="Positive"|Investigation$KTFPPTEPK.classification=="Positive",]
# Investigation <- unique(Investigation)
# Investigation$CTaa[duplicated(Investigation$CTaa)] 
# remove(Investigation)

# Adding the information for which percentage a clone is Dextramer positive.
# Group by clonotype, count cells in that clonotype and cells in that clonotype that have "positive.for.n.dextramers">0 and report percentage of those.
# This has to be done timepoint by timepoint
df_Acute <- as.data.frame(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
df_Recovered <- as.data.frame(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered",])

dt_Acute <- data.table(df_Acute)
dt_Recovered <- data.table(df_Recovered)

dt_Acute <- dt_Acute[,.(dextramer.sum=sum(positive.for.n.dextramers)),by=CTaa] 
dt_Recovered <- dt_Recovered[,.(dextramer.sum=sum(positive.for.n.dextramers)),by=CTaa] 

# Perform first the join of the created data frames
df_Acute <- inner_join(df_Acute,dt_Acute,by = "CTaa")
df_Recovered <- inner_join(df_Recovered,dt_Recovered,by = "CTaa")

rownames(df_Acute) <- df_Acute$Row.names
rownames(df_Recovered) <- df_Recovered$Row.names

df_Acute <- df_Acute[,c(1,48)]
df_Recovered <- df_Recovered[,c(1,48)]

# Perform the second join to the integrated dataset
Integrated_NA_filtered@meta.data <- merge(Integrated_NA_filtered@meta.data,df_Acute,by=0,all.x = TRUE)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names
# colnames(Integrated_NA_filtered@meta.data)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1,3:48,50)]

# head(Integrated_NA_filtered@meta.data)
colnames(df_Recovered) <- c("Row.names.recovered","dextramer.sum.recovered")
Integrated_NA_filtered@meta.data <- merge(Integrated_NA_filtered@meta.data,df_Recovered,by=0,all.x = TRUE)
rownames(Integrated_NA_filtered@meta.data) <- Integrated_NA_filtered@meta.data$Row.names
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1,3:49,51)]
Integrated_NA_filtered@meta.data$dextamer_sum <- coalesce(Integrated_NA_filtered@meta.data$dextramer.sum,Integrated_NA_filtered@meta.data$dextramer.sum.recovered)
Integrated_NA_filtered@meta.data <- Integrated_NA_filtered@meta.data[,c(1:47,50)]

remove(df_Acute,df_Recovered,dt_Acute,dt_Recovered)

Integrated_NA_filtered@meta.data$percent.of.clonotype.positive <- 100*(Integrated_NA_filtered@meta.data$dextamer_sum/Integrated_NA_filtered@meta.data$Frequency)
Integrated_NA_filtered@meta.data$percent.with.values <- paste(100*(Integrated_NA_filtered@meta.data$dextamer_sum/Integrated_NA_filtered@meta.data$Frequency),"(",Integrated_NA_filtered@meta.data$dextamer_sum,"of",Integrated_NA_filtered@meta.data$Frequency,")")

df <- Integrated_NA_filtered@meta.data
only_positives <- df[df$positive.for.n.dextramers>0,]
# colnames(only_positives)
only_positives <- only_positives[,c(1,17,18,19,21,23,24,25,46,49,50)]
write.xlsx(only_positives, 'percentages of clone that is dextramer positive.xlsx')
remove(df,only_positives)


# I enable the option to manually classify clones as Dextramer negative for each single Dextramer and also put the "positive.for.n.dextramers" to 0
# This can be done to filter clones with a chosen "percent.of.clonotype.positive" cutoff.
# Default is no filtering.
# sum(Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0)
# table(Integrated_NA_filtered@meta.data$percent.of.clonotype.positive[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])

Integrated_NA_filtered@meta.data$FTSDYYQLY.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                            Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$TTDPSFLGRY.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                             Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$ATEGALNTPK.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                             Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"
Integrated_NA_filtered@meta.data$KTFPPTEPK.classification[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                            Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- "Negative"

Integrated_NA_filtered@meta.data$positive.for.n.dextramers[Integrated_NA_filtered@meta.data$percent.of.clonotype.positive<0.1 &
                                                             Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0] <- 0

# sum(Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0)
# table(Integrated_NA_filtered@meta.data$percent.of.clonotype.positive[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0],
#       Integrated_NA_filtered@meta.data$Frequency[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])



Integrated@meta.data$Row.names <- rownames(Integrated@meta.data)

# At this point, the preliminary analysis is finished. We use the Seurat objects from this point on to prepare the remaining plots.
# We enable here the option to save and load the  Seurat object for easy loading
SaveH5Seurat(Integrated_NA_filtered, "Integrated_NA_filtered" , overwrite = TRUE)
SaveH5Seurat(Integrated, "Integrated", overwrite = TRUE)


