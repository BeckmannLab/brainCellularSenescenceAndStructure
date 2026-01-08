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


metaFiles=list.files("../meta")
PEC2_sample=fread("PEC2_sample_metadata.txt",data.table=FALSE)
PEC2_mapping=as.data.frame(readxl::read_excel("PEC2_sample_mapping.xlsx"))

#fix for format
PEC2_sample[which(PEC2_sample$Cohort == "Girgenti-snMultiome"),"Cohort"] = "Girgenti-multiome"
PEC2_mapping[which(PEC2_mapping$Cohort == "Girgenti-snMultiome"),"Cohort"] = "Girgenti-multiome"
PEC2_sample[which(PEC2_sample$Individual_ID == "HSB340"),"Individual_ID"] = "HSB340_DFC"
duplicates = c("HSB189", "HSB340", "UMB5297", "RT00389N")

###do for each dataset

#PTSDBrainomics              
folder = "PTSDBrainomics"
expFiles = list.files(paste0(path, "/", folder))
samples = gsub("-annotated_matrix.txt.gz", "", expFiles, fixed = TRUE)
meta_samples = PEC2_sample[PEC2_sample$Cohort == folder, ]
mapping_samples = PEC2_mapping[PEC2_mapping$Cohort == folder, ]

sces = list()
for (expFile in expFiles) {
  sample = gsub("-annotated_matrix.txt.gz", "", expFile, fixed = TRUE)
  if (!sample %in% duplicates) {
    cat("\r", expFile, "\t\t\t\t")
    
    meta_sample = meta_samples[meta_samples$Individual_ID == sample, ]
    mat = fread(paste0(path, "/", folder, "/", expFile), data.table = FALSE)
    rownames(mat) = mat$featurekey
    mat$featurekey = NULL
    
    cellType = colnames(mat)
    meta_sample = meta_sample[rep(seq_len(nrow(meta_sample)), length(cellType)), ]
    meta_sample$cellType = cellType
    row.names(meta_sample) = NULL
    colnames(mat) = paste0(make.names(colnames(mat), unique = TRUE), "_", sample)
    
    # Convert matrix to dgCMatrix instead of spam
    assay_matrix <- as(as.matrix(mat), "dgCMatrix")
    
    sce <- SingleCellExperiment(list(counts = assay_matrix), colData = DataFrame(meta_sample), rowData = DataFrame(feature_id = rownames(mat)))
    sces[[expFile]] = sce
  }
}

sces3 <- do.call(SingleCellExperiment::cbind, sces)
saveRDS(sces3,"sce_all_PTSDBrainomics.RDS")

