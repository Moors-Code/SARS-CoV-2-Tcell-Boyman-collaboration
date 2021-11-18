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
# Part 8 - Creating Donut charts for clonal analysis
###########################################################################################################################################
# The next section creates donut charts showing clone sizes and number of clones for each peptide
# The section can be used to create plots for all peptides, just the peptide need to selected below. Don't forget to change the plot title!
# We are including clones that are found as being Dextramer positive in one disease state, even if they appear as Dex specific only in the other disease state.
Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
Recovered_subset <- subset(Integrated_NA_filtered,subset = Timepoint=="Recovered")
Dex.positive.CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$KTFPPTEPK.classification=="Positive"])
Acute.subset.pos.CTaas <- subset(Acute_subset,subset = CTaa %in% Dex.positive.CTaas) 
Recovered.subset.pos.CTaas <- subset(Recovered_subset, subset = CTaa %in% Dex.positive.CTaas)

Pie_subset_Acute_df <- as.data.frame(table(Acute.subset.pos.CTaas@meta.data$CTaa))
Pie_subset_Recovered_df <- as.data.frame(table(Recovered.subset.pos.CTaas@meta.data$CTaa))

colnames(Pie_subset_Acute_df) <- c("Clonotype","Freq")
colnames(Pie_subset_Recovered_df) <- c("Clonotype","Freq")

Acute_CTaa_df <- Pie_subset_Acute_df
Recovered_CTaa_df <- Pie_subset_Recovered_df

# Now create the donut
# Compute percentages
Acute_CTaa_df$fraction <- Acute_CTaa_df$Freq / sum(Acute_CTaa_df$Freq)
Recovered_CTaa_df$fraction <- Recovered_CTaa_df$Freq / sum(Recovered_CTaa_df$Freq)

# Compute the cumulative percentages (top of each rectangle)
Acute_CTaa_df$ymax = cumsum(Acute_CTaa_df$fraction)
Recovered_CTaa_df$ymax = cumsum(Recovered_CTaa_df$fraction)

# Compute the bottom of each rectangle
Acute_CTaa_df$ymin = c(0, head(Acute_CTaa_df$ymax, n=-1))
Recovered_CTaa_df$ymin = c(0, head(Recovered_CTaa_df$ymax, n=-1))

Acute_CTaa_df$percentage <- Acute_CTaa_df$fraction*100
Recovered_CTaa_df$percentage <- Recovered_CTaa_df$fraction*100

Acute_CTaa_df$rounded <- round(Acute_CTaa_df$percentage,digits = 1)
Recovered_CTaa_df$rounded <- round(Recovered_CTaa_df$percentage,digits = 1)

# Compute label position
Acute_CTaa_df$labelPosition <- (Acute_CTaa_df$ymax + Acute_CTaa_df$ymin) / 2
Recovered_CTaa_df$labelPosition <- (Recovered_CTaa_df$ymax + Recovered_CTaa_df$ymin) / 2

Acute_CTaa_df$Legend <- paste0(Acute_CTaa_df$Clonotype," (",Acute_CTaa_df$rounded,"%)")
Recovered_CTaa_df$Legend <- paste0(Recovered_CTaa_df$Clonotype," (",Recovered_CTaa_df$rounded,"%)")

# Make the plot
ggplot(Acute_CTaa_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Legend)) +
  geom_rect(color='black') +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) + theme_void()+
  geom_text(x = 2, y = 0.5, label = nrow(Acute_CTaa_df), size = 10)+guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle("Acute \n KTFPPTEPK")+
  theme(plot.title = element_text(hjust = 0.5,size=22),legend.position="none")

ggplot(Recovered_CTaa_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Legend)) +
  geom_rect(color='black')  +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) + theme_void()+
  geom_text(x = 2, y = 0.5, label = nrow(Recovered_CTaa_df), size = 10)+guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle("Recovered \n KTFPPTEPK")+
  theme(plot.title = element_text(hjust = 0.5,size=22),legend.position="none")



remove(Acute_CTaa_df, Acute_subset,Acute.subset.pos.CTaas,Pie_subset_Acute,Pie_subset_Acute_df,Pie_subset_Recovered,Pie_subset_Recovered_df,
       Recovered_CTaa_df,Recovered_subset,Recovered.subset.pos.CTaas,Dex.positive.CTaas)

