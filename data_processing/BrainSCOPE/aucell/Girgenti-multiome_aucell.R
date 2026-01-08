ml R/4.4.1
R

library(data.table)
library(SingleCellExperiment)
library(zellkonverter)
library(spam)
library(spam64)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(Seurat)
library(AUCell)
library(dplyr)
library(Matrix)
library(harmony)
library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 48)
library(readxl)
library(GSEABase)
library(HDF5Array)
library(dreamlet)
library(foreach)
library(stringr)
library(tidyr)
library(qs)


# #########
# #AUCell
# #########

sene = readRDS("across_all_celltypes_and_thresholds.RDS") #see across_ds_thresholds
sces3 =readRDS("sce_all_Girgenti-multiome.RDS")  #see set up
sces <- sces3[, colData(sces3)$new_celltype %in% c("Exc", "Oligo", "Int", "Astro", "Micro", "OPC")]
sene_Girgenti_multiome <- sene[sene$Cell %in% colnames(sces), ]

qsave(sene_Girgenti_multiome, "all_Aucell_results_with_cutoffs_Girgenti_multiome.qs")

length(unique(sene_Girgenti_multiome$Cell))
#238101
length(colnames(sces))
#238101

setdiff(unique(sene_Girgenti_multiome$Cell),colnames(sces3))
#0

setdiff(colnames(sces),unique(sene_Girgenti_multiome$Cell))
#0
