CD8+ T cell signature in acute SARS-CoV-2 infection identifies memory precursors

This repository contains the R code that was used to analyze the scRNAseq dataset generated in the "CD8+ T cell signature in acute SARS-CoV-2 infection
identifies memory precursors" project. 
R version: 4.1.0

The file Final_Code_complete.R contains the whole code as one script.
The folder tcell_covid_codes contains the same code divided into 13 parts:
1) preparation
2) pilot_data_addition
3) integration_of_data
4) dimensional_reduction
5) general_condition_comparisons
6) definition_of_specificity
7) covid_specific_cells_plots
8) create_donut_charts
9) persistent_vs_non-persistent_clones
10) GSEA
11) DGEA
12) alluvial_venn_plots
13) additions_to_preprint

Parts 1-6 contain the pre-processing and preparation of the dataset and result in the creation of the Seurat objects "Integrated" and "Integrated_NA_filtered", which are saved at the end of part 6. These parts have to be run one after each other (part 2 requires the output of part 1, ...), where the environment in R should not be cleared between individual parts. 
Parts 7-13 use the two Seurat objects "Integrated" and "Integrated_NA_filtered" to produce most of the scRNAseq analysis plots of the paper. They can be run after loading "Integrated" and "Integrated_NA_filtered".
