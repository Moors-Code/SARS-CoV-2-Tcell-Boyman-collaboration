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
# Part 12 - Alluvial plots, testing if there are clones which are found in more than one SNP cluster, venn diagram
###########################################################################################################################################
# To create alluvial plots, scRepertoire expects as input:
# - a list
# - the list elements should have the names of the samples you want to compare
# I have changed the function "compareClonotypes' of package "Repertoire" so that instead of relative frequencies, absolute frequencies are shown in the alluvial plot.
# trace(compareClonotypes, edit=TRUE)
# Alluvial plot including all persistent clones
persistent_CTaas <- Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Row.names %in% Persistents]
persistent_CTaas <-  unique(persistent_CTaas)

Cells_with_persistent_CTaas <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% persistent_CTaas,]
Acute_Cells_with_persistent_CTaas <- Cells_with_persistent_CTaas[Cells_with_persistent_CTaas$Timepoint=="Acute",]
Recovered_Cells_with_persistent_CTaas <- Cells_with_persistent_CTaas[Cells_with_persistent_CTaas$Timepoint=="Recovered",]

List <- list(Acute_Cells_with_persistent_CTaas,Recovered_Cells_with_persistent_CTaas)
names(List) <- c("Acute","Recovered")

compareClonotypes(List, samples = c("Acute", "Recovered"), 
                  cloneCall="aa", graph = "alluvial")+ggtitle("SARS-CoV-2 common specific clones")+
  theme(plot.title = element_text(hjust = 0.5))

remove(Acute_Cells_with_persistent_CTaas,Cells_with_persistent_CTaas,List,Recovered_Cells_with_persistent_CTaas,persistent_CTaas)



### Testing if there are clones which are found in more than one SNP clusters
# Split by Set1+2 and Set3+4
# head(Integrated_NA_filtered@meta.data)
# table(Integrated_NA_filtered@meta.data$Dataset)
# table(Integrated_NA_filtered@meta.data$assignment)
Cross_patient_clones <- subset(Integrated_NA_filtered, subset = assignment != 'Healthy' & assignment != 'unknown_1' & assignment != 'unknown_2')
Cross_patient_clones_Set1_2 <- subset(Cross_patient_clones, subset = Dataset=="Set1" |Dataset=="Set2")
Cross_patient_clones_Set3_4 <- subset(Cross_patient_clones, subset = Dataset=="Set3" |Dataset=="Set4")

# table(Cross_patient_clones_Set1_2@meta.data$assignment)
# table(Cross_patient_clones_Set3_4@meta.data$assignment)

Cross_patient_clones_Set1_2 <- subset(Cross_patient_clones_Set1_2,subset = positive.for.n.dextramers>0)
# colnames(Cross_patient_clones_Set1_2@meta.data)
Cross_patient_clones@meta.data <- Cross_patient_clones@meta.data[,c(12,19)]
test <- table(Cross_patient_clones@meta.data$CTaa,Cross_patient_clones@meta.data$assignment)
test <- as.data.frame.matrix(test) 
test$zeros <- unname(rowSums(test==0))
# table(test$zeros)
# test[test$zeros==8,]

# Make a pie chart out of this
pie_df <- data.frame(
  Assigned.to.n.SNP.clusters=c("1","2",">2"),
  count=c(4075,159,17))

pie_df <- pie_df %>% 
  arrange(desc(Assigned.to.n.SNP.clusters)) %>%
  mutate(prop = count / sum(pie_df$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie_df$percent <- paste0(pie_df$prop,"%")
pie_df$percent <- c("","95.9%","")
pie_df <- pie_df[order(pie_df$count,decreasing=T),]
pie_df$Assigned.to.n.SNP.clusters <- factor(pie_df$Assigned.to.n.SNP.clusters,levels = c("1", "2", ">2"))
levels(pie_df$Assigned.to.n.SNP.clusters)

ggplot(pie_df, aes(x="", y=prop, fill=Assigned.to.n.SNP.clusters)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(x=1, y = ypos, label = percent), color = "black", size=8) + ggtitle("Number of SNP clusters to which a clone is assigned")+
  theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 20))+ labs(fill = "")

remove(Cross_patient_clones,Cross_patient_clones_Set1_2,Cross_patient_clones_Set3_4,pie_df,test)



### Venn Diagram
Positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0])

Positive_CTaas_in_Acute <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Timepoint=="Acute"&
                                                                          Integrated_NA_filtered@meta.data$CTaa %in% Positive_CTaas])
Positive_CTaas_in_Recovered <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$Timepoint=="Recovered"&
                                                                              Integrated_NA_filtered@meta.data$CTaa %in% Positive_CTaas])
Venn_list <- list('Acute' = Positive_CTaas_in_Acute,'Recovered' = Positive_CTaas_in_Recovered)



# We have noticed that the number of overlapping clones is higher than expected after previous analysis steps. 41 would be expected, but there seem
# to be 44 overlapping clones between the acute and the recovered phase.
# What about the three overlapping clones that should not be there?
# I suspect they are clones that are shared between acute and common but come from different sets! So they are not common clones.
overlap <- intersect(Venn_list[[1]],Venn_list[[2]])
Investigation <- Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% overlap & Integrated_NA_filtered@meta.data$common.clone=="no",]
# unique(Investigation$CTaa) # --> This confirms, that there are 3 non.common clones which add to the 41 "real" overlapping ones.
kicks <- unique(Investigation$CTaa)

# To correct this, I am annotating them in the Positive_CTaas_in_Acute and Positive_CTaas_in_Recovered
Positive_CTaas_in_Acute[grepl(kicks[1],Positive_CTaas_in_Acute)] <- paste0(kicks[1],"_Acute")
Positive_CTaas_in_Acute[grepl(kicks[2],Positive_CTaas_in_Acute)] <- paste0(kicks[2],"_Acute")
Positive_CTaas_in_Acute[grepl(kicks[3],Positive_CTaas_in_Acute)] <- paste0(kicks[3],"_Acute")

Positive_CTaas_in_Recovered[grepl(kicks[1],Positive_CTaas_in_Recovered)] <- paste0(kicks[1],"_Recovered")
Positive_CTaas_in_Recovered[grepl(kicks[2],Positive_CTaas_in_Recovered)] <- paste0(kicks[2],"_Recovered")
Positive_CTaas_in_Recovered[grepl(kicks[3],Positive_CTaas_in_Recovered)] <- paste0(kicks[3],"_Recovered")

Venn_list <- list('Acute' = Positive_CTaas_in_Acute,'Recovered' = Positive_CTaas_in_Recovered)

# Finally, we are creating two different versions of venn digrams, the plot in the figure is a mix of the two.
ggvenn(Venn_list,c('Acute','Recovered'))

Venn_object <- euler(Venn_list)
plot <-plot(Venn_object, alpha=0.6, fill=c("blue", "red"),labels=list(font=1, cex=2))
plot

remove(Investigation, plot,Venn_list,Venn_object,kicks,overlap,Positive_CTaas,Positive_CTaas_in_Acute,Positive_CTaas_in_Recovered)


# Focus on persistent clones
# Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
#                                                                         Integrated_NA_filtered@meta.data$common.clone=="yes"])
# Persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas & 
#                                                            Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
# # How many of these persistents are Dextramer positive
# # table(Integrated_NA_filtered@meta.data$positive.for.n.dextramers[rownames(Integrated_NA_filtered@meta.data)%in%Persistents])
# 
# 
# # These are the cells that belong to persistent clones in the recovered phase
# late_persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas & 
#                                                                 Integrated_NA_filtered@meta.data$Timepoint=="Recovered",])
# # How many of these late_persistents are Dextramer positive
# # table(Integrated_NA_filtered@meta.data$positive.for.n.dextramers[rownames(Integrated_NA_filtered@meta.data)%in%late_persistents])
# 
# remove(late_persistents,Persistents,Common_positive_CTaas)
