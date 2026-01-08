
# -------- DATA DOWNLOAD  --------
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3384nnn/GSM3384109/suppl/GSM3384109_senescence.dge.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3384nnn/GSM3384108/suppl/GSM3384108_HighPDCtrl.dge.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3384nnn/GSM3384107/suppl/GSM3384107_LowPD50Gy.dge.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3384nnn/GSM3384106/suppl/GSM3384106_LowPDCtrl.dge.txt.gz


# -------- LIBRARIES --------
ml R/4.4.1
R

library(data.table)
library(SingleCellExperiment)
library(dreamlet)
library(GSEABase)
library(qvalue)
library(dplyr)
library(assertthat)
library(Matrix)
library(foreach)
library(AUCell)
library(foreach)
library(tidyr)
library(qvalue)
library(stringr)
library("GEOquery")
library(doParallel)
registerDoParallel(cores=10)

# Custom helper functions (edgeR wrappers) - COPY ON MINERVA
source("dreamletCompareClusters_edgeR.R") #see /functions/dreamletCompareClusters_edgeR.R
source("processOneAssay_edgeR.R") #see /functions/processOneAssay_edgeR.R
source("dreamletCompareClusters_edgeR_notPaired.R") #see /functions/dreamletCompareClusters_edgeR_notPaired.R

# -------- FILE PATHS --------

#folders
path_results = "/cutoffs/"
aucell_out = "/aucell/"
expr_results = "/expression/"
de_out = "/DE/"
sen_folder = "/true_senescenceDE/"
true_sen_expression = "/true_senescenceDE/expression/"
true_sen_results = "/true_senescenceDE/DE/"
folder = "/GSE119807/"

#outputs
enrichment_results = "enrichment_or_12.16.25.RDS"
meta_file_save = "meta_for_counts.RDS"
sce_save = "sc_count_with_meta.RDS"
aux_save = "AUCell_GSE119807_1.23.24.RDS"
auc_mat = "all_skin_AUC_mat.RDS"

var_loop = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS"), V2 = c("1percent", 
"1percent", "5percent", "5percent", "10percent", "10percent", 
"20percent", "20percent", "30percent", "30percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_20percent", "pb_20percent", "pb_30percent", 
"pb_30percent"), V4 = c("int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", 
"int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE"
), V5 = c("int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", 
"int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE", 
"exc_FALSE"), V6 = c("int", "exc", "int", "exc", "int", "exc", 
"int", "exc", "int", "exc"), V7 = c(0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4)), class = "data.frame", row.names = c(NA, 
-10L))

# GSE119807_series_matrix.txt  GSM3384106_LowPDCtrl.dge.txt  GSM3384107_LowPD50Gy.dge.txt  GSM3384109_senescence.dge.txt

# Series matrix (sample-level metadata from GEO)
series_matrix = fread(paste0(folder, "GSE119807_series_matrix.txt"))

#Same metadata via GEOquery; exprs(gse) for expression, pData(gse) for pheno
gse=getGEO(filename="/GSE119807_series_matrix.txt")
expression_data <- exprs(gse)
phenotype_data <- pData(gse)

# -------- LOAD SINGLE-CELL COUNT TABLES FOR EACH SAMPLE --------
ctrl1 = fread(paste0(folder, "GSM3384106_LowPDCtrl.dge.txt"))
sen1 = fread(paste0(folder, "GSM3384107_LowPD50Gy.dge.txt"))
sen2 = fread(paste0(folder, "GSM3384109_senescence.dge.txt"))

## add sample-suffixed column names and set gene symbols as rownames
colnames(ctrl1) <- paste0(colnames(ctrl1), "_GSM3384106")
ctrl1 = as.data.frame(ctrl1)
rownames(ctrl1) = ctrl1$GENE_GSM3384106
ctrl1 = ctrl1[,-1]

colnames(sen1) <- paste0(colnames(sen1), "_GSM3384107")
sen1 = as.data.frame(sen1)
rownames(sen1) = sen1$GENE_GSM3384107
sen1 = sen1[,-1]

colnames(sen2) <- paste0(colnames(sen2), "_GSM3384109")
sen2 = as.data.frame(sen2)
rownames(sen2) = sen2$GENE_GSM3384109
sen2 = sen2[,-1]


# -------- MERGE SAMPLES INTO ONE COUNT MATRIX --------
all2 = merge(ctrl1, sen1, by = "row.names")
rownames(all2) = all2$Row.names
all2 = all2[,-1]

all3 = merge(all2, sen2, by = "row.names")
rownames(all3) = all3$Row.names
all3 = all3[,-1]

# Convert to a sparse matrix 
countMat <- as.matrix(all3)
countMat <- as(countMat, "dgCMatrix")

# -------- BUILD COLDATA (SAMPLE LABELS) AND CREATE SCE --------
meta = data.frame(cell = colnames(all3), sample = sub(".*_", "", colnames(all3)), type = "control")
meta$type = "control"
meta$type[meta$sample == "GSM3384109"] = "sen"
meta$type[meta$sample == "GSM3384107"] = "sen"
saveRDS(meta, meta_file_save)
#convert to sce
sce <- SingleCellExperiment(assays = list(counts=countMat), colData=meta)

# saveRDS(sce, sce_save)

# -------- DEFINE GENE SETS (SKIN) AND RUN AUCELL --------
sen_genes = c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
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

skin_genes = c("BRAFV600E","CDKN1A","CDKN2A","CXCL10","CXCL8","CXCR3","H2AJ","H2AX","HMGB1","IL1A","IL1B","IL6","KITLG","LMNB1","MKI67","MMP1","MMP2","MMP3","MMP9","mtDNA4977 del","RARRES2","TNF","TP53BP1")

all = c(sen_genes,skin_genes)
all = unique(all)

#check
length(all)
# [1] 136
assert_that(length(intersect(rownames(all3),all))>0)
#TRUE

#format sen genes
geneSets <- GeneSet(all, setName="sen_list")
geneSets

#aucell
AUCX2 <- AUCell_run(countMat, geneSets=geneSets, BPPARAM=BiocParallel::MulticoreParam(20))

# saveRDS(AUCX2,aux_save)

# -------- THRESHOLD CELLS BY AUCELL SCORE AT MULTIPLE QUANTILES --------
V1 = c("1percent", "5percent", "10percent", "20percent", "30percent")
V2 = c("0.99", "0.95", "0.90", "0.80", "0.70")

all_var =as.data.frame(cbind(V1, V2))
all_var$V2 = as.numeric(all_var$V2)

final_df = c()
resfinal=foreach(i = 1:nrow(all_var)) %do% {
    # subset = subset(meta, major_clust == all_var[i,2])
    # subset_aucx2 = AUCX2[,colnames(AUCX2) %in% rownames(AUCX2) ]
    AUC_mat <- t(getAUC(AUCX2))
    AUC_mat = as.data.frame(AUC_mat)
    AUC_mat$V1 = rownames(AUC_mat)
    threshold <- quantile(AUC_mat$sen_list, all_var[i,2])
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
    AUC_mat$percent = all_var[i,1]
    AUC_mat$cell = rownames(AUC_mat)
    # AUC_mat$sample = sub(".*_", "", AUC_mat$cell)
    # AUC_mat2 = AUC_mat[which(AUC_mat$sample!= "GSM3384108"),] #
    print(i)
    #save
    # saveRDS(AUC_mat, paste0(path_results, all_var[i,2],"_skin_AUC_mat_top_", all_var[i,1],".RDS"))  
    final_df = rbind(AUC_mat, final_df)
}

# saveRDS(final_df, auc_mat)

# -------- PSEUDOBULK AGGREGATION BY (predicted) SENESCENCE LABEL PER SAMPLE --------
all_sene = readRDS(auc_mat)
percentages = c("1percent", "5percent","10percent", "20percent", "30percent")

# Add Sample_ID 
all_sene <-all_sene %>%
  mutate(Sample_ID = str_extract(cell, "(?<=_).*"))
all_sene$type = "skin"
folder = aucell_out

final_df = c()
resfinal=foreach(i = percentages) %do% {
    sene = all_sene[which(all_sene$percent == i),]
    sene$cell_sen = paste0(sene$type,"_",sene$sen)
    sen3 = sene[,c("cell_sen","cell","Sample_ID"), drop = FALSE]
    pb <- list()
    matches <- match(colnames(countMat), sen3$cell)
    valid_matches <- !is.na(matches)
    sce2 <- sce[, valid_matches]
    # colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], , drop = TRUE])
    colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
    colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
    print(i)
    file = paste0(folder, paste0("expression/pseudobulk_", i, ".RDS"))
    if (!file.exists(file)){
        pb = aggregateToPseudoBulk(sce2,
            assay = "counts", 
            cluster_id = "sen_auc",
            sample_id = "Sample_ID",
            BPPARAM = MulticoreParam(30))
        # saveRDS(pb, file)
    }else{
        pb=readRDS(file)
    }
}

# -------- DIFFERENTIAL EXPRESSION (EDGE R VIA DREAMLET WRAPPER) --------

folder = expr_results

var = var_loop
var = var[,1:3]
var = unique(var)


resfinal=foreach(i = 1:nrow(var)) %do% {
      if(length(which(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS")==list.files(de_out,recursive=TRUE)))==0){
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

    ct.pairs <- c("skin_TRUE","skin_FALSE")
    ##run function
    try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, min.samples = 1,min.cells = 1, min.count = 1,method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    #pull
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    # #save
    # saveRDS(fit_top,paste0(de_out,"correct_edger","_",var[i,2], "_fit.RDS"))

}else{
    print(paste(i, "done"))
}}


# -------- SUMMARIZE DE RESULTS ACROSS THRESHOLDS --------
final_df = c()
new_df = c()

for (i in 1:nrow(var)) {
    file_path = paste0(de_out,"correct_edger", "_", var[i,2], "_fit.RDS")
    
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        
        new_df = list()  # Initialize new_df as a list
        new_df$sig = sum(df$FDR <= 0.05)
        new_df$total = nrow(df)
        new_df$percentage = new_df$sig / new_df$total
        new_df$celltype = var[i,6]
        new_df$test = var[i,2]
        new_df$percentage_up = sum(df$logFC > 0) / new_df$total
        new_df$percentage_down = sum(df$logFC < 0) / new_df$total
        new_df$pi1 = 1 - pi0est(df$PValue)$pi0
        
        new_df = as.data.frame(new_df)
        final_df = rbind(final_df, new_df)
    } else {
        cat("File not found:", file_path, "\n")
    }
}


# -------- Enrichment --------

all2 = c()
for (i in 1:nrow(var)){
   file_path = paste0(de_out,"correct_edger", "_", var[i,2], "_fit.RDS")
    
    if (file.exists(file_path)) {
        df = readRDS(file_path)
    # df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}
}
all2$sig = all2$FDR <= 0.05

#sen genes
sen_genes = unique(all)
sen_genes = as.data.frame(unique(sen_genes))
colnames(sen_genes) = "symbol"
dim(sen_genes)
# [1] 136   1

sen_genes$sen_gene = TRUE
intersect(all2$symbol, sen_genes$symbol)
sc_with_sen=merge(all2,sen_genes, by = "symbol", all.x = T)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)]= FALSE

#metrics to loop over
var = var_loop
var2 = var[,c("V3", "V2")]
var = unique(var2)

results <- list()
for (i in 1:nrow(var)) {
    sc_subset <- sc_with_sen[which(sc_with_sen$test == var[i,2]), ]
    contingency_table <- table(sc_subset$sen_gene, sc_subset$sig)
    
    if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
        test <- try(fisher.test(contingency_table), silent = TRUE)
        
        if (!inherits(test, "try-error")) {
            results[[i]] <- list(
                var2 = var[i,2],
                table = contingency_table,
                test = test,
                total_expressed = length(sc_subset$sen_gene),
                degs = sum(sc_subset$FDR <= 0.05),
                kegg_TRUE = sum(sc_subset$sen_gene == TRUE & sc_subset$sig == TRUE),
                total_expressed_kegg = sum(sc_subset$sen_gene == TRUE)
            )
        }
    }
}

final_df <- data.frame(
    var1 = character(),
    var2 = character(),
    p_value = numeric(),
    odds_ratio = numeric(),
    total_expressed = numeric(),
    degs = numeric(),
    kegg_true = numeric(),
    total_expressed_kegg = numeric(),
    stringsAsFactors = FALSE
)

# Iterate through the results list and extract the necessary information
for (i in seq_along(results)) {
    result <- results[[i]]
    
    if (!is.null(result)) {  # Check if the result exists
        # Extract the information and convert to a data frame row
        df_row <- data.frame(
            var1 = var[i, 1],  # Assuming var[i, 1] contains var1
            var2 = result$var2,
            p_value = result$test$p.value,
            odds_ratio = result$test$estimate,
            total_expressed = result$total_expressed,
            degs = result$degs,
            kegg_true = result$kegg_TRUE,
            total_expressed_kegg = result$total_expressed_kegg,
            stringsAsFactors = FALSE
        )
        
        # Bind the data frame row to the final data frame
        final_df <- rbind(final_df, df_row)
    }
}

final_df[order(final_df$var2),]

# -------- fisher test  --------

AUCX2 = readRDS(aux_save)
V1 = c("0.7", "0.8", "0.9", "0.95", "0.99")
V2 = c("30percent", "20percent", "10percent", "5percent", "1percent")
thresholds = as.data.frame(cbind(V1, V2))
thresholds$V1 = as.numeric(thresholds$V1)

df = data.frame(test = character(), odds_ratio = numeric(), fdr = numeric())
for (x in 1:nrow(thresholds)){
    print(x)
    AUC_mat <- t(getAUC(AUCX2))
    AUC_mat = as.data.frame(AUC_mat)
    AUC_mat$V1 = rownames(AUC_mat)
    threshold <- quantile(AUC_mat$sen_list, thresholds[x,1])
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
    AUC_mat$percent = thresholds[i,1]
    AUC_mat$cell = rownames(AUC_mat)

    all = merge(AUC_mat, meta, by.x = "row.names", by.y = "cell")
    all$V4 = 0
    all$V4[all$type == "sen"] = 1
    all$category = "sen"
    all$category[all$type != "sen"] = "control"
    te = fisher.test(table(all$sen, all$category), alternative = "greater")
    df_row <- data.frame(
            test = thresholds[x, 2],  # Assuming var[i, 1] contains var1
            odds_ratio = te$estimate,
            fdr = te$p.value)
    df = rbind(df_row, df)}

# -------- correlate between AUCell scores to true senescence --------

AUC_mat <- t(getAUC(AUCX2))
AUC_mat = as.data.frame(AUC_mat)
AUC_mat$V1 = rownames(AUC_mat)
all = merge(AUC_mat, meta, by.x = "row.names", by.y = "cell")
all$V4 = 0
all$V4[all$type == "sen"] = 1

cor(all$sen_list, all$V4, method = "spearman")
#   0.2922728

# -------- DE correlated to DE they get for that dataset (our method vs there actual defition ) --------

# -------- PSEUDOBULK  --------
sce = readRDS(sce_save)
countMat <-assay(sce)
meta = colData(sce)
meta = as.data.frame(meta)
rownames(meta) = meta$cell
meta = meta[-1,]
meta$sen = FALSE
meta$sen[meta$type == "sen"] = TRUE
meta <- meta %>%
  mutate(Sample_ID = str_extract(cell, "(?<=_).*"))
meta$type = "skin"

folder = sen_folder

sene = meta
sene$cell_sen = paste0(sene$type,"_",sene$sen)
sen3 = sene[,c("cell_sen","cell","Sample_ID"), drop = FALSE]
pb <- list
matches <- match(colnames(countMat), sen3$cell)
valid_matches <- !is.na(matches)
sce2 <- sce[, valid_matches]
colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
file = paste0(folder, paste0("expression/pseudobulk_", "truesen", ".RDS"))
    pb = aggregateToPseudoBulk(sce2,
        assay = "counts", 
        cluster_id = "sen_auc",
        sample_id = "Sample_ID",
        BPPARAM = MulticoreParam(30))
    # saveRDS(pb, file)


# -------- DE with edgeR --------

folder = true_sen_expression

pb = readRDS(paste0(folder, paste0("pseudobulk_", "truesen", ".RDS")))
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

ct.pairs <- c("skin_TRUE","skin_FALSE")
##run function
try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, min.samples = 1,min.cells = 1, min.count = 1,method = "none",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
# pull
fit_top <- as.data.frame(topTags(fit, n = Inf))
print(i)
#save
# saveRDS(fit_top,paste0(true_sen_results,"correct_edger","_","true_sen", "_fit.RDS"))


# -------- pull results --------

df = readRDS(paste0(true_sen_results,"correct_edger","_","true_sen", "_fit.RDS"))
df$symbol = rownames(df)
df$sig = df$FDR <= 0.05
true_sen = df

var = var_loop
var = var[,1:3]
var = unique(var)

all2 = c()
for (i in 1:nrow(var)){
   file_path = paste0(de_out,"correct_edger", "_", var[i,2], "_fit.RDS")
    if (file.exists(file_path)) {
        df = readRDS(file_path)
    # df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}
}

all2$sig = all2$FDR <= 0.05
predicted_sen = all2

##see correlation 
all = merge(true_sen, predicted_sen, by = c("symbol"), all.y = TRUE)

percentages = c("10percent", "20percent", "30percent")
cor_result <- data.frame(test = character(), rho = numeric(), pval = numeric())

for (x in percentages) {
  df <- all[which(all$test == x), ]
  res <- cor.test(df$logFC.x, df$logFC.y, method = "spearman")
  new_row <- data.frame(test = x, rho = res$estimate, pval = res$p.value)
  cor_result <- rbind(cor_result, new_row)
}

cor_result
