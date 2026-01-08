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

#########
#AUCell
#########

sene = readRDS("across_all_celltypes_and_thresholds.RDS") #see across_ds_thresholds
sces3 =readRDS("sce_all_LIBD.RDS") #see set up
sces <- sces3[, colData(sces3)$new_celltype %in% c("Exc", "Oligo", "Int", "Astro", "Micro", "OPC")]
sene_LIBD <- sene[sene$Cell %in% colnames(sces), ]

qsave(sene_LIBD, "all_Aucell_results_with_cutoffs_LIBD.qs")

length(unique(sene_LIBD$Cell))
#49109
length(colnames(sces))
#49109

setdiff(unique(sene_LIBD$Cell),colnames(sces3))
#0

setdiff(colnames(sces),unique(sene_LIBD$Cell))
#0
