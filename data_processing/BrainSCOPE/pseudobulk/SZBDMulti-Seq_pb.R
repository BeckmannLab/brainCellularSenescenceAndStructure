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

#################
#pseudobulk
##############
percentages = c("10percent", "1percent", "20percent", "30percent", "5percent"
    )
final_df = c()
resfinal=foreach(i = percentages) %dopar% {
    for (experiment_name in names(multi_sce)) {
        print(experiment_name)
        name <- str_extract(experiment_name, "^[^-]+")
        sces3 <- qread(paste0(name,"_sce_SZBDMulti-Seq.qs")) #see setup
        all_sene = qread(paste0(name,"_Aucell_results_with_cutoffs_SZBDMulti-Seq.qs")) #see aucell
        all_sene = as.data.frame(all_sene)
        sene = all_sene[which(all_sene$percent == i),]
        sene$cell_sen = paste0(sene$celltype,"_",sene$sen)
        sen3 = sene[,c("cell_sen","Cell","Sample_ID"), drop = FALSE]
        pb <- list()
        matches <- match(colnames(sces3), sen3$Cell)
        valid_matches <- !is.na(matches)
        sce2 <- sces3[, valid_matches]
        colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
        colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
        print(i)
        file = paste0(folder, paste0("SZBDMulti_Seq_pseudobulk_",name,"_", i, ".qs"))
        if (!file.exists(file)){
            pb = aggregateToPseudoBulk(sce2,
                assay = "counts", 
                cluster_id = "sen_auc",
                sample_id = "Sample_ID",
                BPPARAM = MulticoreParam(40))
            qsave(pb, file)
        }else{
            pb=qread(file)
        }
    }
}

