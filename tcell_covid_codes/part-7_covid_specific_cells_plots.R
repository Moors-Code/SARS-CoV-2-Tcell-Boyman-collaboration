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
# Part 7 - Creating plots for the analysis of COVID specific cells
###########################################################################################################################################
# Keep in mind from here on, that there are two Seurat objects: "Integrated" and "Integrated_NA_filtered".
# Most analysis can be done in "Integrated_NA_filtered" alone. 
# But if we want to highlight cells on the UMAP of the full, integrated dataset, we will use "Integrated".

#Integrated_NA_filtered <- LoadH5Seurat("Integrated_NA_filtered.h5Seurat")
#Integrated <- LoadH5Seurat("Integrated.h5Seurat")


# We now progress to highlighting Dextramer positive cells in acute and recovered disease
Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
Recovered_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Recovered")
acute_positives <- rownames(Acute_subset@meta.data[Acute_subset@meta.data$positive.for.n.dextramers>0,])
recovered_positives <- rownames(Recovered_subset@meta.data[Recovered_subset@meta.data$positive.for.n.dextramers>0,])


DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= acute_positives, cols.highlight = 'palevioletred4', cols= "grey92")+
  ggtitle("SARS-CoV-2 peptide specific cells in acute disease")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)+labs(x = "UMAP 1", y="UMAP 2")

DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= recovered_positives, cols.highlight = 'steelblue4', cols= "grey92")+
  ggtitle("SARS-CoV-2 peptide specific cells in recovered disease")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)+labs(x = "UMAP 1", y="UMAP 2")



# As a control we check the patient ID of these cells
# Acute_positives_subset <- subset(Acute_subset,subset = Row.names %in% acute_positives)
# DimPlot(Acute_positives_subset,label = F,group.by = "Patient")
# 
# Recovered_positives_subset <- subset(Recovered_subset,subset = Row.names %in% recovered_positives)
# DimPlot(Recovered_positives_subset,label = F,group.by = "Patient")
# 
# All_positives_subset <- subset(Integrated_NA_filtered, subset = positive.for.n.dextramers>0)
# DimPlot(All_positives_subset,label = F,group.by = "Patient")



remove(Acute_subset,acute_positives,Recovered_subset,recovered_positives)
remove(Acute_positives_subset, All_positives_subset,Recovered_positives_subset)


# A bar plot showing the absolute number of positive cells per cluster and timepoint
df <- Integrated_NA_filtered@meta.data
df <- df[df$Timepoint!="Healthy",]
df <- df %>% group_by(seurat_clusters,Timepoint) %>% add_tally(positive.for.n.dextramers,name = "positive.cells.percluster.and.timepoint") %>% ungroup()

ggplot(df, aes(x=seurat_clusters,y=positive.cells.percluster.and.timepoint, fill=Timepoint))+
  geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0,110))+
  labs(x = "Clusters", y="SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+
  ggtitle("Numbers of Dextramer positive cells per cluster")

remove(df)

# Do the same as percentage plot
df2 <- Integrated_NA_filtered@meta.data
df2 <- df2[df2$Timepoint!="Healthy",]
df2 <- df2 %>% group_by(seurat_clusters) %>% add_count(seurat_clusters,name = "cells.percluster.and.timepoint") %>% ungroup()
df2 <- df2 %>% group_by(seurat_clusters,Timepoint) %>% add_tally(positive.for.n.dextramers,name = "positive.cells.percluster.and.timepoint") %>% ungroup()
df2$percent.pos.cells.per.cluster <- round(100*(df2$positive.cells.percluster.and.timepoint/df2$cells.percluster.and.timepoint),1)

ggplot(df2, aes(x=seurat_clusters,y=percent.pos.cells.per.cluster, fill=Timepoint))+
  geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0,50))+
  labs(x = "Clusters", y="% SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+
  ggtitle("Perentage of Dextramer positive cells per cluster")


remove(df2)

# Create a plot showing in which clusters Dextramer positive cells are in the Acute and Recovered phase.
Contributions <- subset(Integrated_NA_filtered, subset = positive.for.n.dextramers > 0)
Contributions@meta.data$Timepoint <- droplevels(Contributions@meta.data$Timepoint)
cluster.contributions.df <- as.data.frame(with(Contributions@meta.data, table(seurat_clusters, Timepoint)))
cluster.contributions.df <- cluster.contributions.df %>% group_by(Timepoint) %>% add_tally(Freq,name = "total.cells.pertimepoint") %>% ungroup()
cluster.contributions.df <- cluster.contributions.df %>% mutate(percentage.contribution= 100*(Freq/total.cells.pertimepoint))
cluster.contributions.df$percentage.contribution <- round(cluster.contributions.df$percentage.contribution,digits = 1)
cluster.contributions.df$percentage.contribution <- as.numeric(cluster.contributions.df$percentage.contribution)

cluster.contributions.df$Timepoint <- factor(cluster.contributions.df$Timepoint, levels = c("Recovered","Acute","Healthy"))

ggplot(cluster.contributions.df, 
       aes(fill=factor(seurat_clusters, 
                       levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)), 
           y=Timepoint, 
           x=percentage.contribution)) + 
  geom_bar(position="fill", stat="identity",colour="black")+
  theme_classic()+ 
  labs(x = "Cluster contribution", y="Timepoints")+
  scale_x_continuous(labels=scales::percent)+
  labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ 
  theme(text = element_text(size = 17))+
  ggtitle("Distribution of CoV2-Dex+ cells among clusters")+
  ylab("")

remove(cluster.contributions.df,Contributions)
