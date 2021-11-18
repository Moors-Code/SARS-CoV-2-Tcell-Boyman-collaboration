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



###########################################################################################################################################
# Part 10 - GSEA
###########################################################################################################################################
# We want to compare two groups of cells with each other: Persistent vs Nonpersistent cells
##### Input files 
list <- list(Persistents,Non_persistents)
names(list) <- c("persistent", "nonpersistent")

GSEA_persistent.vs.nonpersistent <- subset(Integrated_NA_filtered, 
                                           subset = Row.names %in% list$persistent | Row.names %in% list$nonpersistent)
GSEA_persistent.vs.nonpersistent$persistence <- "nonpersistent"
GSEA_persistent.vs.nonpersistent$persistence[GSEA_persistent.vs.nonpersistent$Row.names %in% list$persistent] <- "persistent"

##### Choose the collection(s) of gene sets to be tested for enrichment 
# see https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# we choose the Hallmark collection:
m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#### Obtain differential expression results with Seurat's FindMarkers 
# We lower the log2 fold change threshold to 0 because
# we want to account also for small changes,
# as GSEA considers all of the genes in an experiment, 
# not only those above an arbitrary cutoff in terms of fold-change or significance.
DefaultAssay(GSEA_persistent.vs.nonpersistent) <- "RNA"
Idents(GSEA_persistent.vs.nonpersistent) <- "persistence"
p_vs_np.results <- FindMarkers(GSEA_persistent.vs.nonpersistent, 
                               logfc.threshold = -Inf,  # equivalent to using 0
                               min.pct = 0.1, 
                               ident.1 = "persistent", ident.2 = "nonpersistent")

#### Define the ranking metric and obtain the pre-ranked feature list 
# For each feature (i.e., each gene) tested, we define the ranking metric as
# the product between the -log of the nominal p-value and the sign of the log2 fold change.
p_vs_np.results <- p_vs_np.results %>% 
  mutate(ranking_metric = -log(p_val)*sign(avg_log2FC) ) %>%
  arrange(desc(ranking_metric)) %>% rownames_to_column(var = "feature")

preranked_list <- p_vs_np.results %>% arrange(desc(ranking_metric)) %>% 
  dplyr::select(feature, ranking_metric) %>% deframe

#### Run GSEA 
# NB: we set the seed before running (F)GSEA as the n permutations (nperm) are random, and we want to make the results reproducible
# The original pre-ranked GSEA (from Subramanian et al., 2005) used 1000 permutations.
# We decide to use more permutations, as FGSEA is fast and allows for that.
# We also decide to run FGSEA-Simple (with permutations),
# rather than FGSEA-multilevel (default option since v3 of the FGSEA pre-print),
# as it's the method consistent with the reference implementation (from Subramanian et al., 2005).
set.seed(42)
fgseaResults <- fgsea(fgsea_sets, stats = preranked_list, nperm=100000) %>%
  as_tibble() %>%
  arrange(desc(NES))

# To print the results to a tab-delimited file, first restructure the list of genes in the leading edge: 
# fgseaResults %>% mutate(leadingEdge = leadingEdge %>% vapply(., paste, collapse = ", ", character(1L))) %>% write.table("gsea_hallmarks_results.txt", sep="\t", quote=F, row.names = F)

#### Plot the NES of the gene sets 
# Only for those found significantly enriched with
# adjusted p-value below a certain threshold (here, padj < 0.1)
p.nes <- ggplot(fgseaResults %>% filter(padj < 0.05), aes(reorder(pathway, NES), NES)) +
  geom_col(fill="blue") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways") + theme_cowplot()

p.nes

#### Save the individual gene set plots 
# Select the gene sets to plot (by p-value adjusted)
genesets_to_plot <- fgseaResults %>% filter(padj < 0.1) %>% select(pathway) %>% unlist %>% as.character

pdf("Fig4f.pdf", onefile=TRUE)
for (i in 1:length(genesets_to_plot)){
  print(plotEnrichment(fgsea_sets[[ genesets_to_plot[i] ]], preranked_list) + 
          # make the gene set names more readable and remove the 'HALLMARK' string at their beginning
          labs(title=stringr::str_replace_all(stringr::str_remove(genesets_to_plot[i], pattern = "HALLMARK_"), 
                                              pattern = "_", replacement = " "), hjust=0.5) +
          theme_cowplot() + labs(x = "Rank", y = "Enrichment Score") +
          geom_line(color = "blue") + geom_point(color = "blue", size=0.1) +
          # center the title
          theme(plot.title = element_text(hjust = 0.5))
  )
}
dev.off()

remove(fgsea_sets,fgseaResults,GSEA_persistent.vs.nonpersistent,list,m_df,p_vs_np.results,p.nes,genesets_to_plot,i,preranked_list)
