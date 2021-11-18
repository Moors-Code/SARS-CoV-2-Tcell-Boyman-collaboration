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
# Part 13 - Additions to preprint version
###########################################################################################################################################
# Fig. 2j
# Comparing cells belonging to persistent clones Acute vs. Recovered (within specific clones!)
Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
                                                                        Integrated_NA_filtered@meta.data$common.clone=="yes"])
dfAcute <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute" & Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas,c(19,21)]
df_Recovered <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered" & Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas,c(19,21)]

colnames(dfAcute)<- c("CTaa","Frequency Acute")
colnames(df_Recovered)<- c("CTaa","Frequency Recovered")

dfAcute <- unique(dfAcute)
df_Recovered <- unique(df_Recovered)

# This dataframe gives an overview over the frequencies that persistent clones have in both timepoints. 
Frequency_trend <- merge(dfAcute,df_Recovered, by="CTaa")
write.xlsx(Frequency_trend, 'Frequencies_of_persistent_clones.xlsx',row.names = TRUE)
persistent_CTaas <- Frequency_trend$CTaa

remove(df_Recovered,dfAcute)

# Get expression levels of a gene for the cells present in a selected group
# Calculate the average expression level
# Make the plot
figure_2_j <- Frequency_trend$CTaa
output2<- matrix(ncol=3, nrow=length(figure_2_j))
i <- 0
while (i<length(figure_2_j)) {
  Acute_expression_levels <- Integrated_NA_filtered@assays$RNA@data["MKI67",rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa==figure_2_j[i+1]&
                                                                                                                        Integrated_NA_filtered@meta.data$Timepoint=="Acute",])]
  Recovered_expression_levels <- Integrated_NA_filtered@assays$RNA@data["MKI67",rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa==figure_2_j[i+1]&
                                                                                                                            Integrated_NA_filtered@meta.data$Timepoint=="Recovered",])]
  Acute_expression_levels <- unname(Acute_expression_levels)
  Recovered_expression_levels <- unname(Recovered_expression_levels)
  Acute_expression_levels <- mean(Acute_expression_levels)
  Recovered_expression_levels <- mean(Recovered_expression_levels)
  
  output2[i+1,1] <- figure_2_j[i+1]
  output2[i+1,2] <- Acute_expression_levels
  output2[i+1,3] <- Recovered_expression_levels
  i <- i+1
}

output2 <- as.data.frame(output2)
colnames(output2) <- c("CTaa","Acute expression level","Recovered expression level")
melted_output2 <- melt(output2,id="CTaa")
melted_output2$value <- as.numeric(melted_output2$value)
melted_output2$value <- round(melted_output2$value,1)
ggplot(melted_output2, aes(x = variable, 
                           y = value,group=CTaa)) +
  geom_line(alpha = 0.8) + 
  geom_point(alpha = 0.7,
             size = 1.5) + 
  theme_classic()


# Creating the excel
write.xlsx(output2, 'MKI67.xlsx',row.names = TRUE)
remove(melted_output2,output2,Acute_expression_levels,i,Recovered_expression_levels,
       figure_2_j)



# Fig. 2h
# Show the localization of individual clones in the UMAP, colored by timepoint
# Generalize this, so that one plot for each persistent clone is generated
i <- 0
while (i<length(persistent_CTaas)) {
  This_round_CTaa <- persistent_CTaas[i+1]
  Acute_Highlight <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute"&
                                                                 Integrated_NA_filtered@meta.data$CTaa %in% This_round_CTaa,])
  Recovered_Highlight <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered"&
                                                                     Integrated_NA_filtered@meta.data$CTaa %in% This_round_CTaa,])
  List <- list(Acute_Highlight,Recovered_Highlight)
  This_rounds_title <- gsub(".*_","",This_round_CTaa)
  names(List) <- c("Acute", "Recovered")
  DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= List, cols.highlight = c("dodgerblue3", "deeppink4"), cols= "grey92")+
    labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 10)+ylim(-10, 15)+ggtitle(This_rounds_title)
  
  pdf(file=paste0("Plot_",This_round_CTaa,".pdf"))
  print(DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= List, cols.highlight = c("dodgerblue3", "deeppink4"), cols= "grey92")+
          labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 10)+ylim(-10, 15)+ggtitle(This_rounds_title))
  
  dev.off()
  i <- i+1
}

remove(List,Acute_Highlight,i,Recovered_Highlight,This_round_CTaa,This_rounds_title)




# Fig. 2i
# Addition on the 14.10., creating our own "gene signature scores" 
active_markers <- list(c("NKG7","PRF1","GZMB","MKI67","CENPU","CENPF"))
resting_markers <- list(c("MT-CO1","MT-CO2","MT-ATP6","IFIT3","IFIT2","TNF"))
Integrated <- AddModuleScore(Integrated,resting_markers,name = "activity.score",assay="RNA")

# Comparing the mean for all cells belonging to persistent clones acute vs recovered
# mean(Integrated@meta.data[Integrated@meta.data$Timepoint=="Acute"&
#                             Integrated@meta.data$CTaa %in% persistent_CTaas,42])
# mean(Integrated@meta.data[Integrated@meta.data$Timepoint=="Recovered"&
#                             Integrated@meta.data$CTaa %in% persistent_CTaas,42])

# Using the AddModuleScore result to make another plot
# For this we need to get rid of the NAs in Integrated
Integrated_noNA <- Integrated
Integrated_noNA@meta.data$CTaa.known <- "yes"
Integrated_noNA@meta.data$CTaa.known[is.na(Integrated_noNA@meta.data$CTaa)] <- "no"
Integrated_noNA <- subset(Integrated_noNA, subset= CTaa.known == "yes")

last_Frequency_trend <- Frequency_trend
figure_2_i <- last_Frequency_trend$CTaa
output2<- matrix(ncol=3, nrow=length(figure_2_i))
i <- 0
while (i<length(figure_2_i)) {
  Acute_score_levels <- Integrated_noNA@meta.data[Integrated_noNA@meta.data$CTaa==figure_2_i[i+1]&
                                                    Integrated_noNA@meta.data$Timepoint=="Acute",42]
  Recovered_score_levels <- Integrated_noNA@meta.data[Integrated_noNA@meta.data$CTaa==figure_2_i[i+1]&
                                                        Integrated_noNA@meta.data$Timepoint=="Recovered",42]
  
  Acute_score_levels <- unname(Acute_score_levels)
  Recovered_score_levels <- unname(Recovered_score_levels)
  Acute_score_levels <- mean(Acute_score_levels)
  Recovered_score_levels <- mean(Recovered_score_levels)
  
  output2[i+1,1] <- figure_2_i[i+1]
  output2[i+1,2] <- Acute_score_levels
  output2[i+1,3] <- Recovered_score_levels
  i <- i+1
}

output2 <- as.data.frame(output2)
colnames(output2) <- c("CTaa","Acute_score_levels","Recovered_score_levels")
melted_output2 <- melt(output2,id="CTaa")
melted_output2$value <- as.numeric(melted_output2$value)
melted_output2$value <- round(melted_output2$value,3)
ggplot(melted_output2, aes(x = variable, 
                           y = value,group=CTaa)) +
  geom_line(alpha = 0.8) + 
  geom_point(alpha = 0.7,
             size = 1.5) + 
  theme_classic()

# Creating the excel
write.xlsx(output2, 'Activity score.xlsx',row.names = TRUE)

remove(active_markers,Integrated_noNA,last_Frequency_trend,melted_output2,output2,resting_markers,Acute_score_levels,i,Recovered_score_levels,figure_2_i)













# Fig. 4g
# Compare expression of selected genes between a persistent and a non-persistent clone
# For the comparison CASSQVIGNQPQHF vs. CASSAPGPLTTQYF
# Creating Excel sheets
genes <- c("B2M","BTG1","CCL4","CENPF","CORO1A","DUSP2","GAPDH","GLRX","HLA-A","HLA-B","HLA-C","IFITM1","MKI67","TXNIP","RRM2","MAP3K8","ZFP36","JUNB","NFKBIA","CENPU","SELL")
gene_counter <- 0
while (gene_counter < length(genes)) {
  
  This_rounds_persistent_cells <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa == "CAANGGDYKLSF_CASSAPGPLTTQYF" &
                                                                              Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
  
  This_rounds_non_persistent_cells <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa == "CALSEHEYNFNKFYF_CASSQVIGNQPQHF" &
                                                                                  Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
  
  This_rounds_persistent_expressions <- Integrated_NA_filtered@assays$RNA@data[genes[gene_counter+1],This_rounds_persistent_cells]
  This_rounds_non_persistent_expressions <- Integrated_NA_filtered@assays$RNA@data[genes[gene_counter+1],This_rounds_non_persistent_cells]
  
  
  df1 <- data.frame(This_rounds_persistent_expressions)
  df1$CTaa <- "persistent_CAANGGDYKLSF_CASSAPGPLTTQYF"
  colnames(df1) <- c("Value","CTaa")
  df2 <- data.frame(This_rounds_non_persistent_expressions)
  df2$CTaa <- "nonpersistent_CALSEHEYNFNKFYF_CASSQVIGNQPQHF"
  colnames(df2) <- c("Value","CTaa")
  df <- rbind(df1,df2)
  
  mean <-aggregate(df[,1], list(df$CTaa), mean)
  colnames(mean)<-c("CTaa", "Value")
  
  write.xlsx(df,paste0(genes[gene_counter+1],'.xlsx'),row.names = TRUE)
  gene_counter <- gene_counter+1  
}

# Which sample are these clones?
Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CAANGGDYKLSF_CASSAPGPLTTQYF",c(11,12,13,25,32:37)]
Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa=="CALSEHEYNFNKFYF_CASSQVIGNQPQHF",c(11,12,13,25,32:37)]

remove(df,df1,df2,mean,gene_counter,genes,This_rounds_non_persistent_cells,This_rounds_non_persistent_expressions,This_rounds_persistent_cells,This_rounds_persistent_expressions)











# Fig. 4f
genes <- c("B2M", "BTG1", "HLA-C", "TXNIP", "CCL4", "JUNB", "RRM2", "CORO1A", "MKI67", "CENPF","CENPU","SELL")
nonpersistent_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & 
                                                                      Integrated_NA_filtered@meta.data$Timepoint=="Acute" &                                                                             Integrated_NA_filtered@meta.data$common.clone=="no"])
genecounter <- 0
output <- matrix(ncol=2, nrow=length(persistent_CTaas))
output2 <- matrix(ncol=2, nrow=length(nonpersistent_CTaas))
while (genecounter < length(genes)) {
  
  CTaa_counter <- 0
  CTaa_counter2 <- 0
  while (CTaa_counter < length(persistent_CTaas)) {
    
    This_round_CTaa <- persistent_CTaas[CTaa_counter+1]
    
    This_rounds_persistent_cells <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa ==This_round_CTaa & Integrated_NA_filtered@meta.data$Timepoint=="Acute",])    
    This_rounds_persistent_expressions<- Integrated_NA_filtered@assays$RNA@data[genes[genecounter+1],This_rounds_persistent_cells]
    mean_expression <- mean(This_rounds_persistent_expressions)
    
    output[CTaa_counter+1,1] <- This_round_CTaa
    output[CTaa_counter+1,2] <- mean_expression
    CTaa_counter <- CTaa_counter+1
  }
  
  while (CTaa_counter2 < length(nonpersistent_CTaas)) {
    
    This_round_CTaa <- nonpersistent_CTaas[CTaa_counter2+1]
    
    This_rounds_nonpersistent_cells <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa ==This_round_CTaa & Integrated_NA_filtered@meta.data$Timepoint=="Acute",])    
    This_rounds_nonpersistent_expressions<- Integrated_NA_filtered@assays$RNA@data[genes[genecounter+1],This_rounds_nonpersistent_cells]
    mean_expression <- mean(This_rounds_nonpersistent_expressions)
    
    output2[CTaa_counter2+1,1] <- This_round_CTaa
    output2[CTaa_counter2+1,2] <- mean_expression
    CTaa_counter2 <- CTaa_counter2+1
  }
  
  persistent_df <- as.data.frame(output)
  colnames(persistent_df) <- c("CTaa", "Value")
  persistent_df$Group <- "persistent"
  nonpersistent_df <- as.data.frame(output2)
  colnames(nonpersistent_df) <- c("CTaa", "Value")
  nonpersistent_df$Group <- "nonpersistent"
  
  
  this_gene_rounds_output <- rbind(persistent_df,nonpersistent_df)
  write.xlsx(this_gene_rounds_output,paste0(genes[genecounter+1],'.xlsx'),row.names = TRUE)
  genecounter <- genecounter+1 
  
}

remove(nonpersistent_df,output,output2,persistent_df,this_gene_rounds_output,CTaa_counter,CTaa_counter2,genecounter,genes,mean_expression,
       This_round_CTaa,This_rounds_nonpersistent_cells,This_rounds_nonpersistent_expressions,This_rounds_persistent_cells,This_rounds_persistent_expressions)




# Ext. Fig. 9
# Differential gene expression and differential protein expression for all cells belonging to persistent clones acute vs recovered (dextr. pos cells only)
Persistent_acutes <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute" & 
                                                                 Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas & 
                                                                 Integrated_NA_filtered@meta.data$positive.for.n.dextramers > 0,])
Persistent_recovereds <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered" & 
                                                                     Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas& 
                                                                     Integrated_NA_filtered@meta.data$positive.for.n.dextramers > 0,])
markers <- as.data.frame(FindMarkers(Integrated_NA_filtered, ident.1 = Persistent_acutes, ident.2 = Persistent_recovereds,assay = "RNA"))
protein.markers <- as.data.frame(FindMarkers(Integrated_NA_filtered, ident.1 = Persistent_acutes, ident.2 = Persistent_recovereds,assay = "Protein"))

remove(Persistent_acutes,Persistent_recovereds,markers,protein.markers)

write.xlsx(markers, 'ext.fig.9.markers.xlsx',row.names = TRUE)
write.xlsx(protein.markers, 'ext.fig.9.proteins.xlsx',row.names = TRUE)


# Fig. 3d
# Comparing cells belonging to persistent clones Acute vs. Recovered (within specific clones!)
Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
                                                                        Integrated_NA_filtered@meta.data$common.clone=="yes"])
dfAcute <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Acute" & Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas,c(19,21)]
df_Recovered <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$Timepoint=="Recovered" & Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas,c(19,21)]

colnames(dfAcute)<- c("CTaa","Frequency Acute")
colnames(df_Recovered)<- c("CTaa","Frequency Recovered")

dfAcute <- unique(dfAcute)
df_Recovered <- unique(df_Recovered)

# This dataframe gives an overview over the frequencies that persistent clones have in both timepoints. 
Frequency_trend <- merge(dfAcute,df_Recovered, by="CTaa")
persistent_CTaas <- Frequency_trend$CTaa

remove(df_Recovered,dfAcute)

# Get expression levels of a gene for the cells present in a selected group
# Calculate the average expression level
# Make the plot
Frequency_trend <- subset(Frequency_trend, Frequency_trend[ , 2] > 2)
Frequency_trend <- subset(Frequency_trend,Frequency_trend_shortened[,3]>2)
figure_3_d <- Frequency_trend$CTaa
output2<- matrix(ncol=3, nrow=length(figure_3_d))
i <- 0
while (i<length(figure_3_d)) {
  Acute_expression_levels <- Integrated_NA_filtered@assays$Protein@data["CCR7.1",rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa==figure_3_d[i+1]&
                                                                                                                             Integrated_NA_filtered@meta.data$Timepoint=="Acute",])]
  Recovered_expression_levels <- Integrated_NA_filtered@assays$Protein@data["CCR7.1",rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa==figure_3_d[i+1]&
                                                                                                                                 Integrated_NA_filtered@meta.data$Timepoint=="Recovered",])]
  Acute_expression_levels <- unname(Acute_expression_levels)
  Recovered_expression_levels <- unname(Recovered_expression_levels)
  Acute_expression_levels <- mean(Acute_expression_levels)
  Recovered_expression_levels <- mean(Recovered_expression_levels)
  
  output2[i+1,1] <- figure_3_d[i+1]
  output2[i+1,2] <- Acute_expression_levels
  output2[i+1,3] <- Recovered_expression_levels
  i <- i+1
}

output2 <- as.data.frame(output2)
colnames(output2) <- c("CTaa","Acute expression level","Recovered expression level")
melted_output2 <- melt(output2,id="CTaa")
melted_output2$value <- as.numeric(melted_output2$value)
melted_output2$value <- round(melted_output2$value,1)
ggplot(melted_output2, aes(x = variable, 
                           y = value,group=CTaa)) +
  geom_line(alpha = 0.8) + 
  geom_point(alpha = 0.7,
             size = 1.5) + 
  theme_classic()


# Creating the excel
write.xlsx(output2, 'CCR7.1.xlsx',row.names = TRUE)
remove(melted_output2,output2,Acute_expression_levels,i,Recovered_expression_levels,
       figure_3_d)

remove(melted_output2,output2,Acute_expression_levels,i,Recovered_expression_levels,figure_3_d)


