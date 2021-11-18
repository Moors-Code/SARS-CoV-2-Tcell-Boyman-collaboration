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
# Part 5 - General analysis comparing cells from Healthy, Acute and Recovered patients and Severity conditions
###########################################################################################################################################
Integrated$Timepoint <- factor(x = Integrated$Timepoint, levels = c("Healthy", "Acute","Recovered"))
coul2 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2","maroon", "deeppink1", "blue1", "steelblue4","green1","yellow3")
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)

# We noticed that clusters 12 and 13 are very small and not interesting for the analysis, therefore we exclude them at this point
Integrated <- subset(Integrated, idents = c(0:11))
Integrated@meta.data$seurat_clusters <- droplevels(Integrated@meta.data$seurat_clusters)
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.1,group.by = "Patient")+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)

# Creating a dot plot
features <- c("HLA-DRB5","HLA-DPB1","HLA-DQB1","EOMES","TIGIT","RRM2","TK1","CENPF","CENPM","MKI67","MCM4","CENPU",
              "IFIT1","IFIT2","IFIT3","IFIH1","TNF","LTA","IFNG","GZMB","GZMH","GNLY","PRF1","NKG7",
              "TCF7", "SELL", "IL7R", "CCR7")
DotPlot(Integrated, features = features,assay = "RNA",cols = c("lightgrey", "orangered1"))+ xlab("Marker genes")+ylab("Cluster")+coord_flip()
remove(features)

# Creating a UMAP where cells are colored by SNP cluster
# Prepare this also as percentage barplot
# table(Integrated@meta.data$assignment)
Integrated_w.o_unkown <- subset(Integrated,subset = assignment !="unknown_1" & assignment !="unknown_2" & assignment !="Healthy")
# table(Integrated_w.o_unkown@meta.data$assignment)
DimPlot(Integrated_w.o_unkown, reduction = "umap",label = F,pt.size = 0.1,group.by = "assignment")+labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("")

df <- Integrated_w.o_unkown@meta.data
df <- df %>% group_by(seurat_clusters,assignment) %>% add_count(seurat_clusters,name = "cells.percluster.and.assignment") %>% ungroup()
df <- df %>% group_by(assignment) %>% add_count(assignment,name = "cells.per.assignment") %>% ungroup()
df$percent.cells.per.cluster <- round(100*(df$cells.percluster.and.assignment/df$cells.per.assignment),1)
# colnames(df)
df <- df[,c(12,39,42)]
df <- unique(df)

# Every row in this bar plot is one of 20 SNP clusters (the two timepoints are not separated).
ggplot(df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11)), y=assignment, x=percent.cells.per.cluster)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "", y="SNP assignment")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(df,Integrated_w.o_unkown)


# Now I am creating individual UMAPs for the three timepoints, only one color per plot
Healthy <- subset(Integrated, subset = Timepoint=="Healthy")
Acute <- subset(Integrated, subset = Timepoint=="Acute")
Recovered <- subset(Integrated, subset = Timepoint=="Recovered")

DimPlot(Healthy, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'gray40')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Healthy")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)

DimPlot(Acute, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'palevioletred4')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Acute")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)

DimPlot(Recovered, reduction = "umap",group.by = 'orig.ident',label = F,pt.size = 0.25, cols = 'steelblue4')+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Recovered")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 10)+ylim(-10, 15)


# Show the cluster contributions per timepoint as a stacked percentage barplot
cluster.contributions.df <- as.data.frame(with(Integrated@meta.data, table(seurat_clusters, Timepoint)))
cluster.contributions.df <- cluster.contributions.df %>% group_by(Timepoint) %>% add_tally(Freq,name = "total.cells.pertimepoint") %>% ungroup()
cluster.contributions.df <- cluster.contributions.df %>% mutate(percentage.contribution= 100*(Freq/total.cells.pertimepoint))
cluster.contributions.df$percentage.contribution <- round(cluster.contributions.df$percentage.contribution,digits = 1)
cluster.contributions.df$percentage.contribution <- as.numeric(cluster.contributions.df$percentage.contribution)

cluster.contributions.df$Timepoint <- factor(cluster.contributions.df$Timepoint, levels = c("Recovered","Acute","Healthy"))

ggplot(cluster.contributions.df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)), y=Timepoint, x=percentage.contribution)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Cluster contribution", y="Timepoints")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(cluster.contributions.df)
remove(Acute, Healthy, Recovered)

# Add severity as column to the object
Integrated@meta.data$Severity <- "unkown"
mild_patients <- c('Patient 21','Patient 25', 'Patient 27','Patient 35','Patient 38','Patient 120','Patient 133')
severe_patients <- c('Patient 1', 'Patient 5', 'Patient 33', 'Patient 41', 'Patient 45', 'Patient 48', 'Patient 50',
                     'Patient 57','Patient 63','Patient 68','Patient 103','Patient 106','Patient 112')
Integrated@meta.data$Severity[Integrated@meta.data$Patient%in%mild_patients] <- "mild"
Integrated@meta.data$Severity[Integrated@meta.data$Patient%in%severe_patients] <- "severe"

remove(mild_patients,severe_patients)


# Individual UMAPs for the acute severities
Mild <- subset(Integrated, subset = Severity=="mild")
Severe <- subset(Integrated, subset = Severity=="severe")
Mild_acute <- subset(Mild, subset = Timepoint=="Acute")
Severe_acute <- subset(Severe, subset = Timepoint=="Acute")
Mild_recovered <- subset(Mild, subset = Timepoint=="Recovered")
Severe_recovered <- subset(Severe, subset = Timepoint=="Recovered")

DimPlot(Mild_acute, reduction = "umap",group.by = 'seurat_clusters',label = F,pt.size = 0.1,cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Mild acute illness")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 15)+ylim(-10, 15)

DimPlot(Severe_acute, reduction = "umap",group.by =  'seurat_clusters',label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Severe acute illness")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 15)+ylim(-10, 15)

DimPlot(Mild_recovered, reduction = "umap",group.by = 'seurat_clusters',label = F,pt.size = 0.1,cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Mild_recovered illness")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 15)+ylim(-10, 15)

DimPlot(Severe_recovered, reduction = "umap",group.by =  'seurat_clusters',label = F,pt.size = 0.1, cols = coul2)+labs(x = "UMAP 1", y="UMAP 2")+
  ggtitle("Severe_recovered illness")+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.position="none")+xlim(-10, 15)+ylim(-10, 15)

remove(Mild,Mild_acute,Mild_recovered,Severe,Severe_acute,Severe_recovered)

