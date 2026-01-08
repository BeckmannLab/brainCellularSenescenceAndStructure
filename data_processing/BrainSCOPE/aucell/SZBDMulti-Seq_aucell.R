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
sces =readRDS("sce_all_SZBDMulti-Seq.RDS") #see set up
sene <- sene %>%
      mutate(Sample_ID = str_extract(Cell, "(?<=_).*"))
multi_sce = sces 

# Iterate over each SingleCellExperiment object
for (experiment_name in names(multi_sce)) {
  sces3 <- multi_sce[[experiment_name]]
  name <- str_extract(experiment_name, "^[^-]+")
  
  cat("Processing experiment:", experiment_name, "\n")
  
  colData(sces3)$new_celltype <- ifelse(colData(sces3)$cellType %in% 
                                      c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", 
                                        "L5/6 NP", "L6 IT", "L6 IT Car3", 
                                        "L6 CT", "L6b"), "Exc", 
                                      ifelse(colData(sces3)$cellType %in% 
                                             c("Pvalb", "Sst", "Sst Chodl", "Lamp5", 
                                               "Lamp5 Lhx6", "Chandelier", "Vip", 
                                               "Sncg", "Pax6"), "Int", 
                                             colData(sces3)$cellType))

    sces3 <- sces3[, colData(sces3)$new_celltype %in% c("Exc", "Oligo", "Int", "Astro", "Micro", "OPC")]

    sene_SZBDMulti_Seq <- sene[sene$Cell %in% colnames(sces3), ]
    assert_that(length(unique(sene_SZBDMulti_Seq$Cell)) == length(colnames(sces3)))
    qsave(sene_SZBDMulti_Seq, paste0(name,"_Aucell_results_with_cutoffs_SZBDMulti-Seq.qs"))
    qsave(sces3, paste0(name,"_sce_SZBDMulti-Seq.qs"))
}

