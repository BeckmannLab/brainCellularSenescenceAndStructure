
ml R/4.4.1
R

rm(list=ls())
options(stringsAsFactors=F)
library(data.table)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(Seurat)
library(AUCell)
library(dplyr)
library(Matrix)
library(harmony)
library(BiocParallel)
library(readxl)
Sys.setenv(OMP_NUM_THREADS = 48)
library(GSEABase)
library(HDF5Array)
library(dreamlet)
library(doParallel)
registerDoParallel(cores=40)
library(foreach)
library(stringr)
library(tidyr)

###run DE on AU_cell

setwd("path")

lbp = #LBP single cell RNAseq data

#subset to liv cohort
ids <- c("PT-0236R", "PT-0237R", "PT-0256L", "PT-0257R", "PT-0261L", 
"PT-0262L", "PT-0264L", "PT-0264R", "PT-0267L", "PT-0267R", "PT-0280L", 
"PT-0280R", "PT-0282L", "PT-0282R", "PT-0283L", "PT-0284R", "PT-0285L", 
"PT-0285R", "PT-0286L", "PT-0286R", "PT-0287L", "PT-0288L", "PT-0289L", 
"PT-0289R", "PT-0290L", "PT-0290R", "PT-0292L", "PT-0292R", "PT-0295L", 
"PT-0296L", "PT-0297L")

celltype = c( "Oli", "Int", "Exc", "Ast", "MG", "OPC")

##ids
tnbc2 <- subset(lbp, subset = Sample_ID %in% ids)
liv_data =GetAssayData(tnbc2, slot="counts", assay="RNA")

#senescence genes
genes_42 <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", 
"CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", 
"GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", 
"LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", 
"MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", 
"STING1", "TGFB1", "TIMP2", "TNF", "TP53")
genes_125 <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
"BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", 
"CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", 
"CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", 
"CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", 
"CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", 
"FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", 
"HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", 
"IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", 
"IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", 
"INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", 
"MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", 
"NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", 
"PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", 
"SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", 
"TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", 
"TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")

all = c(genes_42,genes_125)
all = unique(all)

geneSets <- GeneSet(all, setName="sen_list")
geneSets

var = structure(list(V1 = c("10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent", "1percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent", 
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "20percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "30percent", "30percent"
), V2 = c("MG", "Exc", "Oli", "Int", "NonNeu", "Ast", "OPC", 
"MG", "Exc", "Oli", "Int", "NonNeu", "Ast", "OPC", "MG", "Exc", 
"Oli", "Int", "NonNeu", "Ast", "OPC", "MG", "Exc", "Oli", "Int", 
"NonNeu", "Ast", "OPC", "MG", "Exc", "Oli", "Int", "NonNeu", 
"Ast", "OPC"), V3 = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.99, 
0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.95, 0.95, 0.95, 0.95, 0.95, 
0.95, 0.95, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 
0.7, 0.7, 0.7, 0.7)), row.names = c(NA, -35L), class = "data.frame")

#run AUCell
AUCX2 <- AUCell_run(liv_data, 
                     geneSets=geneSets,
                     BPPARAM=BiocParallel::MulticoreParam(20))

final_df = c()
resfinal=foreach(i = 1:nrow(var)) %do% {
    subset = subset(tnbc2, CellTypeAWC == var[i,2])
    subset_aucx2 = AUCX2[,colnames(AUCX2) %in% colnames(subset) ]
    AUC_mat <- t(getAUC(subset_aucx2))
    AUC_mat = as.data.frame(AUC_mat)
    threshold <- quantile(AUC_mat$sen_list, var[i,3])
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
    AUC_mat$celltype = var[i,2]
    AUC_mat$percent = var[i,1]
    AUC_mat$cell = rownames(AUC_mat)
    print(i)
    #save
    saveRDS(AUC_mat, paste0("/path/", var[i,2],"_146_AUC_mat_top_", var[i,1],".RDS"))  
    final_df = rbind(AUC_mat, final_df)
}

#pull results
all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("/path/", var[i,2],"_146_AUC_mat_top_", var[i,1],".RDS"))
    all2 = rbind(df,all2)
}

saveRDS(all2, "all_Aucell_results_with_cutoffs.RDS")

############
#pseudobulk
#############
all_sene = readRDS("all_Aucell_results_with_cutoffs.RDS")

percentages = c("10percent", "1percent", "20percent", "30percent", "5percent"
    )
all_sene <- all_sene %>%
        mutate(Sample_ID = str_extract(cell, "^[^_]+"))

folder = "path"
sce <- as.SingleCellExperiment(lbp)

final_df = c()
resfinal=foreach(i = percentages) %do% {
    sene = all_sene[which(all_sene$percent == i),]
    sene$cell_sen = paste0(sene$celltype,"_",sene$sen)
    sen3 = sene[,c("cell_sen","cell","Sample_ID"), drop = FALSE]
    pb <- list()
    matches <- match(colnames(sce), sen3$cell)
    valid_matches <- !is.na(matches)
    sce2 <- sce[, valid_matches]
    colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], , drop = TRUE])
    colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
    colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
    print(i)
    file = paste0(folder, paste0("expression2/pseudobulk_", i, ".RDS"))
    if (!file.exists(file)){
        pb = aggregateToPseudoBulk(sce2,
            assay = "counts", 
            cluster_id = "sen_auc",
            sample_id = "Sample_ID",
            BPPARAM = MulticoreParam(40))
        saveRDS(pb, file)
    }else{
        pb=readRDS(file)
    }
}


####run DE
folder = "path"

var = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS"
), V2 = c("1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent"
), V4 = c("Ast_TRUE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", 
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE", "Ast_FALSE", 
"Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", 
"OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_TRUE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", 
"NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE"), V5 = c("Ast_FALSE", "Exc_FALSE", 
"Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", 
"Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", 
"Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE", 
"MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", "Ast_TRUE", 
"Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", 
"OPC_FALSE", "Ast_FALSE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", 
"NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE"), V6 = c("ast", "exc", 
"int", "mg", "noneu", "oli", "opc", "ast", "exc", "int", "mg", 
"noneu", "oli", "opc", "ast", "exc", "int", "mg", "noneu", "oli", 
"opc", "ast", "exc", "int", "mg", "noneu", "oli", "opc", "ast", 
"exc", "int", "mg", "noneu", "oli", "opc"), V7 = c(0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)), row.names = c(NA, -35L), class = "data.frame")

resfinal=foreach(i = 1:nrow(var)) %do% {
      if(length(which(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS")==list.files("path",recursive=TRUE)))==0){
    pb = readRDS(paste0(folder,var[i,1]))
    ncells = as.data.frame(int_colData(pb)$n_cells)
    ncells$celltype = rownames(ncells) 

    ncells2 = ncells %>%
      pivot_longer(!celltype,names_to = "id", values_to = "count")

    ncells2 = as.data.frame(ncells2)

    ncells2$id = gsub("\\.", "_", ncells2$id)

    colnames(ncells2) = c("sen_auc", "Sample_ID", "ncell")

    ##add ncells

    metadata(pb)$aggr_means$Indv_ID = metadata(pb)$aggr_means$Sample_ID 

    metadata(pb)$aggr_means$Indv_ID=gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)

    metadata(pb)$aggr_means$ncells = as.numeric(ncells2$ncell)
    # sce_list[[var$name[i]]] = pb

    ct.pairs <- c(var[i,4], var[i,5])
    ##run function
    try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    #pull
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    #save
    saveRDS(fit_top,paste0("path/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))

}else{
    print(paste(i, "done"))
}}


################
#save results
###############

library(limma)
final_df = c()
new_df = c()

for (i in 1:nrow(var)){
    df = readRDS(paste0("path/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    new_df$sig = length(which(df$FDR <= 0.05))
    sig = length(which(df$FDR <= 0.05))
    new_df$total = nrow(df)
    total = nrow(df)
    new_df$percentage = (sig / total) 
    new_df$celltype = var[i,6]
    new_df$test = var[i,2]
    length_up = sum(df$logFC > 0)
    length_down = sum(df$logFC < 0)
    new_df$percentage_up = (length_up/ total)
    new_df$percentage_down = (length_down/ total)
    new_df$pi1=1 - propTrueNull(df$PValue)
    new_df = as.data.frame(new_df)
    final_df = rbind(final_df, new_df)
    }

saveRDS(final_df, "DE_results_across.RDS")
