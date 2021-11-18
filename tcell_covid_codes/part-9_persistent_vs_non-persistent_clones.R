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
# Part 9 - Comparing persistent vs non-persistent clones
###########################################################################################################################################
# The first plot creates a pie chart that shows all acute clones and whether they are also found in the 6 month timepoint or not
# Acute_subset <- subset(Integrated_NA_filtered, subset = Timepoint=="Acute")
# unique_acute_CTaas <- unique(Acute_subset@meta.data$CTaa)
# Acute_commons <- Acute_subset@meta.data[Acute_subset@meta.data$common.clone=="yes",]
# unique_acute_commons <- unique(Acute_commons$CTaa)
# 
# persistent_df <- data.frame(
#   Occurance=c("persistent","non-persistent"),
#   count=c(length(unique_acute_commons),(length(unique_acute_CTaas)-length(unique_acute_commons))))
# 
# persistent_df <- persistent_df %>% 
#   arrange(desc(Occurance)) %>%
#   mutate(prop = count / sum(persistent_df$count) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )
# 
# ggplot(persistent_df, aes(x="", y=prop, fill=Occurance)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) +
#   theme_void() +
#   geom_text(aes(y = ypos, label = count), color = "white", size=6) +
#   scale_fill_manual(values = c("dodgerblue3", "deeppink4")) +
#   ggtitle("Persistent vs. non-persistent clones \n (All unique clones)") +
#   theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 26))
# 
# remove(Acute_commons,Acute_subset,persistent_df,unique_acute_commons,unique_acute_CTaas)


# Create the same plot only for COVID specific clones (clones that are specific in acute phase but not necessarily specific in recovered phase)!
# New, include also those that are found to be positive in the recovered phase but not necessarily positive in the acute phase. 
Acute_but_non_common_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & 
                                                                             Integrated_NA_filtered@meta.data$Timepoint=="Acute" &
                                                                             Integrated_NA_filtered@meta.data$common.clone=="no"]) 
Common_positive_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
                                                                        Integrated_NA_filtered@meta.data$common.clone=="yes"])

persistent_df2 <- data.frame(
  Occurance=c("persistent","non-persistent"),
  count=c(length(Common_positive_CTaas),length(Acute_but_non_common_CTaas)))

persistent_df2 <- persistent_df2 %>% 
  arrange(desc(Occurance)) %>%
  mutate(prop = count / sum(persistent_df2$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(persistent_df2, aes(x="", y=prop, fill=Occurance,colo)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(y = ypos, label = count), color = "white", size=12) +
  scale_fill_manual(values = c("dodgerblue3", "deeppink4")) +ggtitle("Persistent vs. non-persistent clones \n (Specific unique clones in Acute phase)")+
  theme(plot.title = element_text(hjust = 0.5,size=22),text = element_text(size = 26))

remove(persistent_df2)


# Highlight specific persistent and non persistent clones in the UMAP

###
#best_specifics <- intersect(unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & Integrated_NA_filtered@meta.data$Timepoint=="Acute"]),
#                                              unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 & Integrated_NA_filtered@meta.data$Timepoint=="Recovered"])
#)
#`%notin%` <- Negate(`%in%`)
#Non_persistents_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
#                                                                  Integrated_NA_filtered@meta.data$CTaa %notin% best_specifics])
#Persistents_CTaas <- unique(Integrated_NA_filtered@meta.data$CTaa[Integrated_NA_filtered@meta.data$positive.for.n.dextramers>0 &
#                                                                    Integrated_NA_filtered@meta.data$CTaa %in% best_specifics])                            
#
#
#Non_persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Non_persistents_CTaas&
#                                                               Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
#Persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Persistents_CTaas&
#                                                               Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
#
###

Non_persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Acute_but_non_common_CTaas& 
                                                               Integrated_NA_filtered@meta.data$Timepoint=="Acute",])
Persistents <- rownames(Integrated_NA_filtered@meta.data[Integrated_NA_filtered@meta.data$CTaa %in% Common_positive_CTaas & 
                                                           Integrated_NA_filtered@meta.data$Timepoint=="Acute",])


list <- list(Persistents,Non_persistents)
names(list) <- c("Persistent", "Non-persistent")
DimPlot(Integrated, label=F,repel = T, group.by="seurat_clusters", cells.highlight= list, cols.highlight = c("deeppink4", "dodgerblue3"), cols= "grey92")+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 10)+ylim(-10, 15)+ggtitle("")

remove(list,Acute_but_non_common_CTaas,Common_positive_CTaas)


# To which clusters do the persistent and non persistent clones belong?
ps_and_nps <- subset(Integrated_NA_filtered,subset = Row.names %in% Persistents |Row.names %in% Non_persistents)
ps_and_nps@meta.data <- ps_and_nps@meta.data %>% group_by(common.clone,seurat_clusters) %>% add_count(seurat_clusters,name = "cells.per.cluster") %>% ungroup()

# sum(ps_and_nps@meta.data$common.clone=="yes")
# sum(ps_and_nps@meta.data$common.clone=="no")

ps_and_nps@meta.data$percent <- 0
ps_and_nps@meta.data$percent[ps_and_nps@meta.data$common.clone =="yes"] <- round(100*(ps_and_nps@meta.data$cells.per.cluster[ps_and_nps@meta.data$common.clone =="yes"]/sum(ps_and_nps@meta.data$common.clone=="yes")),1)
ps_and_nps@meta.data$percent[ps_and_nps@meta.data$common.clone =="no"] <- round(100*(ps_and_nps@meta.data$cells.per.cluster[ps_and_nps@meta.data$common.clone =="no"]/sum(ps_and_nps@meta.data$common.clone=="no")),1)

df <- ps_and_nps@meta.data
# colnames(df)
df <- df[,c(23,39,52)]
df <- unique(df)
ggplot(df, aes(x=seurat_clusters,y=percent, fill=common.clone))+geom_bar(position="dodge", stat="identity",colour="black")+theme_classic()+
  scale_fill_manual(values=c("palevioletred4", "steelblue4"))+ scale_y_continuous(expand = c(0, 0), limits = c(0,50))+labs(x = "Clusters", y="% SARS-CoV-2 peptide specific cells")+
  theme(text = element_text(size = 15))+labs(fill='Persistent clone')+ggtitle("Distribution of cells belonging to Persistent \n or Non-persistent clones across clusters")

ggplot(df, aes(fill=factor(seurat_clusters, levels=c(0,1,2,3,4,5,6,7,8,9,10,11)), y=common.clone, x=percent)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "", y="Persistent")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cluster")+
  scale_fill_manual(values=c(coul2))+ theme(text = element_text(size = 17))

remove(df,ps_and_nps)


# Violin plot comparing clone sizes specific vs non-specific
Violin.plot.persistents <- Integrated_NA_filtered@meta.data[Persistents,]
Violin.plot.nonpersistents <- Integrated_NA_filtered@meta.data[Non_persistents,]

# colnames(Violin.plot.persistents)
# colnames(Violin.plot.nonpersistents)
Violin.plot.persistents$occurence <- "Persistents"
Violin.plot.nonpersistents$occurence <- "Non-persistents"

Violin.plot.df <- rbind(Violin.plot.persistents,Violin.plot.nonpersistents)
# colnames(Violin.plot.df)
Violin.plot.df <- Violin.plot.df[,c(19,21,51)]
Violin.plot.df <- unique(Violin.plot.df)

ggplot(Violin.plot.df, aes(x=occurence, y=Frequency, fill=occurence)) + 
  geom_violin()+theme_classic()+geom_jitter(shape=16, position=position_jitter(0.0))+
  scale_fill_manual(values=c("dodgerblue3", "deeppink4"))+ labs(x="",y="Clonal Frequency")+ theme(legend.title = element_blank(),text = element_text(size = 20))+
  scale_y_continuous(trans='log2')

remove(Violin.plot.df,Violin.plot.nonpersistents,Violin.plot.persistents)



# Highlighting the 10 most expanded persistent and nonpersistent clones in the UMAP -> This has to be adapted based on the cutoff chosen in Part 6.
# First persistent
Specifics_persistent <- Integrated_NA_filtered
Specifics_persistent@meta.data$Row.names <- rownames(Specifics_persistent@meta.data)
Specifics_persistent@meta.data$Specific <- "no"
Specifics_persistent@meta.data$Specific[rownames(Specifics_persistent@meta.data)%in%Persistents] <- "yes" # this has to be changed depending on the group
Specifics_persistent@meta.data$Clone <- "a"
Specifics_persistent@meta.data$Clone[Specifics_persistent@meta.data$Specific=="yes"] <- Specifics_persistent@meta.data$CTaa[Specifics_persistent@meta.data$Specific=="yes"]
Specifics_persistent <- subset(Specifics_persistent, subset = Specific =="yes")
Specifics_persistent@meta.data$Clone <- as.numeric(as.factor(Specifics_persistent@meta.data$Clone))
Specifics_persistent@meta.data <- Specifics_persistent@meta.data %>% group_by(CTaa) %>% add_count(Clone, name = "fake.Clonesize") %>% ungroup()
Specifics_persistent@meta.data <- as.data.frame(Specifics_persistent@meta.data)
rownames(Specifics_persistent@meta.data) <- Specifics_persistent@meta.data$Row.names
# table(Specifics_persistent@meta.data$fake.Clonesize)

First <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==83, ])
Second <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==77, ])
Third <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==40, ])
Fourth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==27, ])
Fifth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==26, ])
Sixth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==22, ])
Seventh <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==16, ])
Eigth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==12, ])
Nineth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==10, ])
Tenth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==9, ])

Integrated@meta.data$Highlight.column <- "AAA"
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%First] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%First])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Second] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Second])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Third] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Third])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fourth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fourth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fifth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fifth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Sixth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Sixth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Seventh] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Seventh])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Eigth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Eigth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Nineth] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Nineth])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Tenth][c(1,3,6,7,8,9,10,12,15)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Tenth])[1]

highlight.colors <- c("grey92",coul2[1:10])
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.8, group.by = "Highlight.column",order = T,cols = highlight.colors)+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("Top 10 expanded persistent clones")+theme(text = element_text(size = 10))


# Then the non-persistent
Specifics_persistent <- Integrated_NA_filtered
Specifics_persistent@meta.data$Row.names <- rownames(Specifics_persistent@meta.data)
Specifics_persistent@meta.data$Specific <- "no"
Specifics_persistent@meta.data$Specific[rownames(Specifics_persistent@meta.data)%in%Non_persistents] <- "yes" # this has to be changed depending on the group
Specifics_persistent@meta.data$Clone <- "a"
Specifics_persistent@meta.data$Clone[Specifics_persistent@meta.data$Specific=="yes"] <- Specifics_persistent@meta.data$CTaa[Specifics_persistent@meta.data$Specific=="yes"]
Specifics_persistent <- subset(Specifics_persistent, subset = Specific =="yes")

Specifics_persistent@meta.data$Clone <- as.numeric(as.factor(Specifics_persistent@meta.data$Clone))
Specifics_persistent@meta.data <- Specifics_persistent@meta.data %>% group_by(CTaa) %>% add_count(Clone, name = "fake.Clonesize") %>% ungroup()
# table(Specifics_persistent@meta.data$fake.Clonesize)
Specifics_persistent@meta.data <- as.data.frame(Specifics_persistent@meta.data)
rownames(Specifics_persistent@meta.data) <- Specifics_persistent@meta.data$Row.names
# table(Specifics_persistent@meta.data$fake.Clonesize)

First <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==11, ])
Second <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==9, ])
Third <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==5, ])
Fourth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==5, ])
Fifth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Sixth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Seventh <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Eigth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==4, ])
Nineth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==3, ])
Tenth <- rownames(Specifics_persistent@meta.data[Specifics_persistent@meta.data$fake.Clonesize==3, ])

Integrated@meta.data$Highlight.column <- "AAA"
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%First] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%First])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Second] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Second])
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Third][c(1,2,3,4,5)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Third])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fourth][c(6,7,8,9,10)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fourth])[2]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Fifth][c(1,2,4,7)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Fifth])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Sixth][c(3,6,10,13)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Sixth])[2]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Seventh][c(5,8,9,15)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Seventh])[3]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Eigth][c(11,12,14,16)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Eigth])[4]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Nineth][c(1,2,3)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Nineth])[1]
Integrated@meta.data$Highlight.column[rownames(Integrated@meta.data)%in%Tenth][c(4,5,6)] <- unique(Integrated@meta.data$CTaa[rownames(Integrated@meta.data)%in%Tenth])[2]

highlight.colors <- c("grey92",coul2[1:10])
DimPlot(Integrated, reduction = "umap",label = F,pt.size = 0.8, group.by = "Highlight.column",order = T,cols = highlight.colors)+
  labs(x = "UMAP 1", y="UMAP 2")+xlim(-10, 15)+ylim(-10, 15)+ggtitle("Top 10 expanded non-persistent clones")+theme(text = element_text(size = 10))

remove(Specifics_persistent,First,Second,Third,Fourth,Fifth,Sixth,Seventh,Eigth,Nineth,Tenth,highlight.colors,coul2)
