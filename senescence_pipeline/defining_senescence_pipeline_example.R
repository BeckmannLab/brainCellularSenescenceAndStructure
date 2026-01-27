
# ————————————————————————————————————————————————————————————————————————————
# 0. Setup
# ————————————————————————————————————————————————————————————————————————————

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = 48)

# Load packages, one per line
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
library(GSEABase)
library(HDF5Array)
library(dreamlet)
library(foreach)
library(doParallel)
library(stringr)
library(tidyr)
library(limma)

# Base directory
setwd()

#matrix
lbp = #seurat object with all counts, etc. 

# ————————————————————————————————————————————————————————————————————————————
# 1. Subset to samples of interest and format (dataset specific)
# ————————————————————————————————————————————————————————————————————————————

sample_ids <- c("S00021","S00253","S00324","S00556","S00579","S00935","S01447","S01462",
  "S11088","S17633","S17680","S19335","T-111","T-164","T-168","T-2011",
  "T-4639.3","T-4915","T-5354","T-563")
tnbc2     <- subset(lbp, subset = Sample_ID %in% sample_ids)
pm_counts <- GetAssayData(tnbc2, slot = "counts", assay = "RNA")

# ————————————————————————————————————————————————————————————————————————————
# 2. AUCell scoring
# ————————————————————————————————————————————————————————————————————————————

# 2.1 Define gene sets - brain specific (modify based on organ)
genes_42  <- c("4-HNE","AXL","BCL2","CCL2","CCL3","CCL4","CCL5","CDKN1A",
"CDKN2A","CDKN2B","CDKN2D","CSF1","CSF2RA","CXCL1","CXCL8",
"GLB1","H2AX","HMGB1","IGF1","IL1A","IL1B","IL27","IL6",
"LGALS3","LGALS3BP","LMNB1","MACROH2A1","MIF","MMP12","MMP3",
"MTOR","PCNA","PLAUR","SA-β-Gal","SATB1","SERPINE1","SPP1",
"STING1","TGFB1","TIMP2","TNF","TP53") #Saul et al. 2022

genes_125 <- c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3",
"BMP2","BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20",
"CCL24","CCL26","CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8",
"CD55","CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB",
"CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8",
"CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2",
"FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF",
"HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3",
"IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15",
"IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7",
"INHA","IQGAP2","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF",
"MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3","MMP9",
"NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU",
"PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG",
"SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX",
"TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B",
"TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2") #Suryadevara et al. 2024

all_genes <- unique(c(genes_42, genes_125))
geneSets  <- GeneSetCollection(GeneSet(all_genes, setName = "sen_list"))

# 2.2 Run AUCell
AUCX2 <- AUCell_run(pm_counts, geneSets = geneSets,
                    BPPARAM = MulticoreParam(20))

# 2.3 Threshold & save per celltype/percent

#variables
var = readRDS("../var_aucell_final.RDS")

#loop
final_df <- foreach(i = seq_len(nrow(var)), 
                    .combine = rbind,      # row-bind each AUC_mat
                    .packages = c("stringr")) %do% {

  # pull parameters for this iteration
  pct       <- var[i, 1]
  cell_type <- var[i, 2]
  thr_prop  <- var[i, 3]

  # subset cells & AUC object
  cells_sub   <- subset(tnbc2, CellTypeAWC == cell_type)
  auc_sub     <- AUCX2[, colnames(AUCX2) %in% colnames(cells_sub)]

  # compute AUC matrix
  AUC_mat <- t(getAUC(auc_sub)) %>% 
    as.data.frame()

  # threshold and annotate
  cutoff <- quantile(AUC_mat$sen_list, thr_prop)
  AUC_mat$sen       <- AUC_mat$sen_list >= cutoff
  AUC_mat$celltype  <- cell_type
  AUC_mat$percent   <- pct
  AUC_mat$cell      <- rownames(AUC_mat)

  # save per-celltype file
  fname <- paste0(cell_type, "_146_AUC_mat_top_", pct, ".RDS")
  # saveRDS(AUC_mat, file.path("../au_cell/pm", fname))

  message("done iteration ", i, ": ", cell_type, " @ ", pct)

  # return it for rbind
  AUC_mat
}

# saveRDS(final_df, "../all_Aucell_results_with_cutoffs.RDS")

# ————————————————————————————————————————————————————————————————————————————
# 3. Pseudobulk 
# ————————————————————————————————————————————————————————————————————————————

sce     <- as.SingleCellExperiment(lbp)
all_sene <- final_df %>%
  mutate(Sample_ID = str_extract(cell, "^[^_]+"))

percentages <- c("10percent","1percent","20percent","30percent","5percent")

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
        # saveRDS(pb, file)
    }else{
        pb=readRDS(file)
    }
}

# ————————————————————————————————————————————————————————————————————————————
# 4. Differential expression (edgeR via dreamlet)
# ————————————————————————————————————————————————————————————————————————————
registerDoParallel(cores = 10)

#modified dreamletCompareClusters - should be open access
source("/sc/arion/projects/psychgen/lbp/data/neuroimaging/anina/au_cell/functions/dreamletCompareClusters_edgeR.R")
source("/sc/arion/projects/psychgen/lbp/data/neuroimaging/anina/au_cell/functions/processOneAssay_edgeR.R")


#var to loop over
var <- readRDS("../var_DE_final.RDS")

#4.1 DE loop for each cut off and cell type
resfinal=foreach(i = 1:nrow(var)) %do% {
      if(length(which(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS")==list.files("../DE/",recursive=TRUE)))==0){
    #load
    pb = readRDS(paste0(folder,var[i,1]))
    #format
    ncells = as.data.frame(int_colData(pb)$n_cells)
    ncells$celltype = rownames(ncells) 
    ncells2 = ncells %>%
      pivot_longer(!celltype,names_to = "id", values_to = "count")
    ncells2 = as.data.frame(ncells2)
    ncells2$id = gsub("\\.", "_", ncells2$id)
    colnames(ncells2) = c("sen_auc", "Sample_ID", "ncell")
    metadata(pb)$aggr_means$Indv_ID = metadata(pb)$aggr_means$Sample_ID 
    metadata(pb)$aggr_means$Indv_ID=gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
    metadata(pb)$aggr_means$ncells = as.numeric(ncells2$ncell)
    #define
    ct.pairs <- c(var[i,4], var[i,5])
    ##run function
    try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    #pull
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    #save
    # saveRDS(fit_top,paste0("../DE/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))

}else{
    print(paste(i, "done"))
}}

# ————————————————————————————————————————————————————————————————————————————
# 5. Summarize & save overall DE metrics (extra step not necessary)
# ————————————————————————————————————————————————————————————————————————————

final_df = c()
new_df = c()

for (i in 1:nrow(var)){
    df = readRDS(paste0("../DE/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    # if (length(which(df$FDR <= 0.05)==0){
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

# saveRDS(final_df, "../DE_results_across.RDS")

# ————————————————————————————————————————————————————————————————————————————
# 6. Chosing cutoffs based on enrichment
# ————————————————————————————————————————————————————————————————————————————

#6.1 pull all results
all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("../DE/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

all2$sig = all2$FDR <= 0.05

#format gene list
all_genes = as.data.frame(unique(all_genes))
colnames(all_genes) = "symbol"
dim(all_genes)
# [1] 146   1

all_genes$sen_gene = TRUE

#check
intersect(all2$symbol, all_genes$symbol)

sc_with_sen=merge(all2,all_genes, by = "symbol", all.x = T)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)]= FALSE

#6.2 enrichments
final_df <- do.call(rbind, lapply(seq_len(nrow(var)), function(i) {
  #select
  pct       <- var[i, "V2"]
  cell_type <- var[i, "V6"]   
  #subset
  sc_sub <- sc_with_sen[
    sc_with_sen$test     == pct &
    sc_with_sen$celltype == cell_type, 
  ]
  
  ct <- with(sc_sub, table(sen_gene, sig))
  if (any(dim(ct) < 2)) return(NULL)
  
  #Fisher’s test
  ft <- try(fisher.test(ct), silent=TRUE)
  if (inherits(ft, "try-error")) return(NULL)
  
  #output
  data.frame(
    percent               = pct,
    celltype              = cell_type,
    p_value               = ft$p.value,
    odds_ratio            = unname(ft$estimate),
    total_expressed       = nrow(sc_sub),
    degs                  = sum(sc_sub$FDR <= 0.05),
    kegg_true             = sum(sc_sub$sen_gene & sc_sub$sig),
    total_expressed_kegg  = sum(sc_sub$sen_gene),
    stringsAsFactors      = FALSE
  )
}))

#sort
final_df <- final_df[order(final_df$celltype), ]

####For each cell type, you should:

# 1. Compare all cutoffs by their odds‐ratio.
# 2. Pick the cutoff that has the highest odds‐ratio.
# 3. If two (or more) cutoffs tie on odds‐ratio, choose the one with the smallest p-value.
# 4. Proceed with downstream analyses using the DE results from that cell-type & cutoff combination.For example, by these rules ast ends up using the 10% cutoff, while exc uses the 1% cutoff.

# > final_df
#      percent celltype      p_value  odds_ratio total_expressed  degs kegg_true
# 1   1percent      ast 5.134584e-13   23.128680           12495   177        13
# 8   5percent      ast 3.275242e-29   34.229622           10602   733        34
# 15 20percent      ast 1.693493e-29   30.637493           11507  1137        39
# 22 30percent      ast 4.231148e-31   36.329062           11278  1319        42
# 29 10percent      ast 2.420762e-31   41.617111           10364   878        37
# 2   1percent      exc 2.554625e-16   23.520010           13955  5949        52
# 9   5percent      exc 8.121261e-09   11.687286           14216  8375        50
# 16 20percent      exc 1.259584e-06    6.212562           15600 10051        56
# 23 30percent      exc 1.914020e-06    5.505289           16291 10263        56
# 30 10percent      exc 6.042019e-08   10.194943           14771  9571        56
