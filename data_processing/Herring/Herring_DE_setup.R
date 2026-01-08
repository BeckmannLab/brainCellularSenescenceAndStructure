ml R/4.4.1
R

# load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(dreamlet)
library(zenith)
library(DelayedArray)
library(GSEABase)
library(muscat)
library(cowplot)
library(scater)
library(tidyverse)
library(kableExtra)
library(qvalue)
library(scattermore)
library(RColorBrewer)
library(corrplot)
library(viridis)
library(dplyr)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(ggrepel)
library(Seurat)
library(data.table)
library(assertthat)
library(AUCell)
library(matrixStats)
library(doParallel)
library(foreach)
registerDoParallel(cores=40)


sce = readH5AD("Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad", use_hdf5=TRUE, verbose=TRUE, reader="R") #see methods for download info

meta=fread("Processed_data_RNA-all_BCs-meta-data.csv",data.table=FALSE) #see methods for download info

geneMeta=fread("Processed_data_RNA-all_genes-meta-data.csv",data.table=FALSE) #see methods for download info

#format
rownames(meta)=meta$V1
meta$V1=NULL
colData(sce)=DataFrame(meta)
rownames(geneMeta)=geneMeta$V1
geneMeta$V1=NULL
rowData(sce)=DataFrame(geneMeta)
assay(sce,1) <- as(assay(sce,1), "dgCMatrix")
countMat=assay(sce, 1)
cpmMat=assay(sce, 2)

#checks
summary(meta$doublet_score)
summary(meta$percent_mito)
summary(meta$percent_ribo)
table(meta$cell_type)
table(meta$major_clust)
table(meta$sub_clust)

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

#check
assert_that(length(intersect(rownames(geneMeta),all))>0)

geneSets <- GeneSet(all, setName="sen_list")

AUCX2 <- AUCell_run(countMat, geneSets=geneSets, BPPARAM=BiocParallel::MulticoreParam(20))

meta2 <- meta %>%
  mutate(major_clust_summarized = case_when(
    major_clust %in% c("PN_dev","L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4") ~ "exc",
    major_clust %in% c("MGE_dev", "CGE_dev", "ID2", "VIP", "SST", "PV","LAMP5_NOS1", "PV_SCUBE3") ~ "int",
    TRUE ~ major_clust 
  ))


all_var = structure(list(V1 = c("1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"10percent", "10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "30percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent"), V2 = c("Astro", "CGE_dev", "ID2", "L2-3_CUX2", 
"L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5_NOS1", "MGE_dev", 
"Micro", "Oligo", "OPC", "PN_dev", "PV", "Vas", "exc", "int", 
"Astro", "CGE_dev", "ID2", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", 
"L5-6_TLE4", "LAMP5_NOS1", "MGE_dev", "Micro", "Oligo", "OPC", 
"PN_dev", "PV", "Vas", "exc", "int", "Astro", "CGE_dev", "ID2", 
"L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5_NOS1", 
"MGE_dev", "Micro", "Oligo", "OPC", "PN_dev", "PV", "Vas", "exc", 
"int", "Astro", "CGE_dev", "ID2", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", 
"L5-6_TLE4", "LAMP5_NOS1", "MGE_dev", "Micro", "Oligo", "OPC", 
"PN_dev", "PV", "Vas", "exc", "int", "Astro", "CGE_dev", "ID2", 
"L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5_NOS1", 
"MGE_dev", "Micro", "Oligo", "OPC", "PN_dev", "PV", "Vas", "exc", 
"int"), V3 = c(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 
0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.95, 0.95, 
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 
0.95, 0.95, 0.95, 0.95, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 
0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 
0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 
0.7, 0.7, 0.7, 0.7)), row.names = c(NA, -85L), class = "data.frame")


#define thresholds
final_df = c()
resfinal=foreach(i = 1:nrow(all_var)) %do% {
    subset = subset(meta, major_clust == all_var[i,2])
    subset_aucx2 = AUCX2[,colnames(AUCX2) %in% rownames(subset) ]
    AUC_mat <- t(getAUC(subset_aucx2))
    AUC_mat = as.data.frame(AUC_mat)
    threshold <- quantile(AUC_mat$sen_list, as.numeric(all_var[i,3]))
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
    AUC_mat$celltype = all_var[i,2]
    AUC_mat$percent = all_var[i,1]
    AUC_mat$cell = rownames(AUC_mat)
    print(i)
    #save
    saveRDS(AUC_mat, paste0(all_var[i,2],"_146_AUC_mat_top_", all_var[i,1],".RDS"))  
    final_df = rbind(AUC_mat, final_df)
}

#pull together
all2 = c()
for (i in 1:nrow(all_var)){
    df = readRDS(paste0(all_var[i,2],"_146_AUC_mat_top_", all_var[i,1],".RDS"))
    all2 = rbind(df,all2)
}
saveRDS(all2, "all_aucell_cutoffs_8.1.24.RDS")

###do for grouping
all_var = structure(list(V1 = c("1percent", "1percent", "5percent", "5percent", 
"10percent", "10percent", "20percent", "20percent", "30percent", 
"30percent"), V2 = c("exc", "int", "exc", "int", "exc", "int", 
"exc", "int", "exc", "int"), V3 = c(0.99, 0.99, 0.95, 0.95, 0.9, 
0.9, 0.8, 0.8, 0.7, 0.7)), row.names = c(NA, -10L), class = "data.frame")

final_df = c()
resfinal=foreach(i = 1:nrow(all_var)) %do% {
    subset = subset(meta2, major_clust_summarized == all_var[i,2])
    subset_aucx2 = AUCX2[,colnames(AUCX2) %in% rownames(subset) ]
    AUC_mat <- t(getAUC(subset_aucx2))
    AUC_mat = as.data.frame(AUC_mat)
    threshold <- quantile(AUC_mat$sen_list, as.numeric(all_var[i,3]))
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
    AUC_mat$celltype = all_var[i,2]
    AUC_mat$percent = all_var[i,1]
    AUC_mat$cell = rownames(AUC_mat)
    print(i)
    #save
    saveRDS(AUC_mat, paste0(all_var[i,2],"_146_AUC_mat_top_", all_var[i,1],".RDS"))  
    final_df = rbind(AUC_mat, final_df)
}

#########
#pseudobulk
#############

all_sene = readRDS("all_aucell_cutoffs_8.1.24.RDS")
percentages = c("10percent", "1percent", "20percent", "30percent", "5percent")
all_sene <- all_sene %>%
        mutate(Sample_ID = str_extract(cell, "(?<=-).*"))

#pseudobulk for each cell type and cut off
final_df = c()
resfinal=foreach(i = percentages) %do% {
    sene = all_sene[which(all_sene$percent == i),]
    sene$cell_sen = paste0(sene$celltype,"_",sene$sen)
    sen3 = sene[,c("cell_sen","cell","Sample_ID"), drop = FALSE]
    pb <- list()
    matches <- match(colnames(countMat), sen3$cell)
    valid_matches <- !is.na(matches)
    sce2 <- sce[, valid_matches]
    colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
    colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
    print(i)
    file = paste0(folder, paste0("expression/pseudobulk_", i, ".RDS"))
    if (!file.exists(file)){
        pb = aggregateToPseudoBulk(sce2,
            assay = "X", 
            cluster_id = "sen_auc",
            sample_id = "Sample_ID",
            BPPARAM = MulticoreParam(30))
        saveRDS(pb, file)
    }else{
        pb=readRDS(file)
    }
}
