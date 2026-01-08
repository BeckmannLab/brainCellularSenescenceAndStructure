# -------- DATA DOWNLOAD --------
# wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175533/suppl/GSE175533%5Fsceasy%5Fhay%2Eh5ad%2Egz"

# -------- LIBRARIES --------
options(stringsAsFactors = FALSE)

library(zellkonverter)          # readH5AD
library(SingleCellExperiment)
library(dreamlet)
library(GSEABase)
library(qvalue)
library(dplyr)
library(Matrix)
library(AUCell)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(cores = 10)

# dreamlet/edgeR helpers 
source("dreamletCompareClusters_edgeR.R") #see functions
source("processOneAssay_edgeR.R") #see functions
source("dreamletCompareClusters_edgeR_notPaired.R") #see functions

# -------- FILE PATHS --------
base                <- "/GSE175533/"
aucell_out          <- file.path(base, "aucell/")
cutoff_dir          <- file.path(aucell_out, "cutoffs/")
expr_dir            <- file.path(aucell_out, "expression/")
de_dir              <- file.path(aucell_out, "DE/")
true_sen_root       <- file.path(aucell_out, "true_senescenceDE/")
true_expr_dir       <- file.path(true_sen_root, "expression/")
true_de_dir         <- file.path(true_sen_root, "DE/")

# outputs
meta_rds            <- file.path(base, "meta_Data.RDS")
aux_rds             <- file.path(base, "AUCell_GSE175533_1.16.24.RDS")
all_auc_mat_rds     <- file.path(cutoff_dir, "all_lung_AUC_mat.RDS")
enrich_rds          <- file.path(aucell_out, "enrichment_or_12.16.25.RDS")

# ensure directories exist
invisible(lapply(c(aucell_out, cutoff_dir, expr_dir, de_dir, true_sen_root, true_expr_dir, true_de_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# -------- LOAD H5AD & PICK ASSAY SAFELY --------
sce <- readH5AD(file.path(base, "GSE175533_sceasy_hay.h5ad"))
# choose an assay: prefer "counts", else "X", else first available
assay_name <- if ("counts" %in% assayNames(sce)) "counts" else
              if ("X" %in% assayNames(sce)) "X" else assayNames(sce)[1]
countMat   <- assay(sce, assay_name)

# meta frame with stable columns
meta <- as.data.frame(colData(sce))
meta$barcode <- rownames(meta)     # explicit barcode column
# this dataset uses 'exp' with values like "WT" and "tert"
stopifnot("exp" %in% colnames(meta))
# saveRDS(meta, meta_rds)

# -------- DEFINE GENE SETS (LUNG) --------
sen_genes <- c(
  "ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3",
  "CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL3L1","CCL4",
  "CCL5","CCL7","CCL8","CD55","CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1",
  "CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","CXCR2",
  "DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7",
  "GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2",
  "IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18",
  "IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2",
  "ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14",
  "MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT",
  "PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F",
  "SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C",
  "TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"
)
lung_genes <- c("BCL2L1","CDKN1A","CDKN2A","CXCL8","GDF15","H2AX","IL10","IL1A",
                "IL1B","IL6","MMP12","MMP8","SERPINE1","TGFB1","TNF","TNFRSF1B","TP53","VEGFA")
all_genes  <- unique(c(sen_genes, lung_genes))
stopifnot(length(intersect(rownames(countMat), all_genes)) > 0)
geneSets   <- GeneSet(all_genes, setName = "sen_list")

# -------- RUN AUCELL --------
AUCX2 <- AUCell_run(countMat, geneSets = geneSets, BPPARAM = BiocParallel::MulticoreParam(20))
# saveRDS(AUCX2, aux_rds)

# -------- THRESHOLD CELLS BY AUCELL SCORE --------
percent_labels <- c("1percent","5percent","10percent","20percent","30percent")
quantiles      <- c(0.99, 0.95, 0.90, 0.80, 0.70)
all_var        <- data.frame(percent = percent_labels, q = quantiles)

auc_list <- vector("list", nrow(all_var))
for (i in seq_len(nrow(all_var))) {
  AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
  AUC_mat$cell <- rownames(AUC_mat)
  thr <- quantile(AUC_mat$sen_list, all_var$q[i])
  AUC_mat$sen <- AUC_mat$sen_list >= thr
  AUC_mat$percent <- all_var$percent[i]
  q_lbl <- gsub("\\.", "_", sprintf("%.2f", all_var$q[i]))
  # saveRDS(AUC_mat, file.path(cutoff_dir, sprintf("%s_lung_AUC_mat_top_%s.RDS", q_lbl, all_var$percent[i])))
  auc_list[[i]] <- AUC_mat
}
all_sene <- do.call(rbind, auc_list)
# saveRDS(all_sene, all_auc_mat_rds)

# -------- PSEUDOBULK (predicted labels) --------
# Sample_ID: try last token after "-" (common), else fall back to 'exp' (WT/tert)
meta$Sample_ID <- str_extract(meta$barcode, "(?<=-)[^-_]+$")
meta$Sample_ID[is.na(meta$Sample_ID) | meta$Sample_ID == ""] <- as.character(meta$exp)
meta$type <- "lung"

percentages <- percent_labels
pb_var <- data.frame(
  file  = paste0("pseudobulk_", percentages, ".RDS"),
  test  = percentages,
  pb    = paste0("pb_", percentages),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pb_var))) {
  sene <- subset(all_sene, percent == pb_var$test[i])
  sene$cell_sen <- paste0("lung_", sene$sen)
  sen3 <- sene[, c("cell_sen","cell"), drop = FALSE]

  matches <- match(colnames(countMat), sen3$cell)
  valid   <- !is.na(matches)
  sce2    <- sce[, valid]
  colData(sce2)$sen_auc   <- as.character(sen3[matches[valid], "cell_sen", drop = TRUE])
  # bring Sample_ID from meta by column order of sce2
  colData(sce2)$Sample_ID <- meta$Sample_ID[match(colnames(sce2), meta$barcode)]

  out_file <- file.path(expr_dir, pb_var$file[i])
  if (!file.exists(out_file)) {
    pb <- aggregateToPseudoBulk(
      sce2, assay = assay_name,
      cluster_id = "sen_auc",
      sample_id  = "Sample_ID",
      BPPARAM = BiocParallel::MulticoreParam(30)
    )
    # saveRDS(pb, out_file)
  }
}

# -------- DIFFERENTIAL EXPRESSION (edgeR via dreamlet wrapper) --------
for (i in seq_len(nrow(pb_var))) {
  in_file  <- file.path(expr_dir, pb_var$file[i])
  out_file <- file.path(de_dir, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(in_file) || file.exists(out_file)) next

  pb <- readRDS(in_file)

  # add ncells metadata
  ncells  <- as.data.frame(int_colData(pb)$n_cells)
  ncells$celltype <- rownames(ncells)
  ncells2 <- ncells %>% pivot_longer(!celltype, names_to = "id", values_to = "count") %>% as.data.frame()
  ncells2$id <- gsub("\\.", "_", ncells2$id)
  colnames(ncells2) <- c("sen_auc", "Sample_ID", "ncell")

  metadata(pb)$aggr_means$Indv_ID <- metadata(pb)$aggr_means$Sample_ID
  metadata(pb)$aggr_means$Indv_ID <- gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
  metadata(pb)$aggr_means$ncells  <- as.numeric(ncells2$ncell)

  ct.pairs <- c("lung_TRUE","lung_FALSE")
  fit <- try(dreamletCompareClusters_edgeR(
    pb, ct.pairs,
    min.samples = 1, min.cells = 1, min.count = 1,
    method = "fixed", min.prop = 0.4,
    useProcessOneAssay_edgeR = TRUE
  ), silent = TRUE)

  if (!inherits(fit, "try-error")) {
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    # saveRDS(fit_top, out_file)
  }
}

# -------- SUMMARIZE DE RESULTS --------
sum_list <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_dir, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(f)) next
  df <- readRDS(f)
  sum_list[[i]] <- data.frame(
    sig              = sum(df$FDR <= 0.05),
    total            = nrow(df),
    percentage       = sum(df$FDR <= 0.05) / nrow(df),
    celltype         = "lung",
    test             = pb_var$test[i],
    percentage_up    = sum(df$logFC > 0) / nrow(df),
    percentage_down  = sum(df$logFC < 0) / nrow(df)
  )
}
summary_df <- do.call(rbind, sum_list)
print(summary_df[order(match(summary_df$test, percentages)), ])

# -------- Enrichment of senescence genes among DE --------
all_de <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_dir, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(f)) next
  df <- readRDS(f)
  df$test   <- pb_var$test[i]
  df$symbol <- rownames(df)
  df$sig    <- df$FDR <= 0.05
  all_de[[i]] <- df
}
all_de <- do.call(rbind, all_de)

sen_tbl <- data.frame(symbol = unique(all_genes), sen_gene = TRUE)
sc_with_sen <- merge(all_de, sen_tbl, by = "symbol", all.x = TRUE)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)] <- FALSE

enrich_list <- list()
for (i in seq_len(nrow(pb_var))) {
  sc_subset <- subset(sc_with_sen, test == pb_var$test[i])
  tab <- table(sc_subset$sen_gene, sc_subset$sig)
  if (nrow(tab) >= 2 && ncol(tab) >= 2) {
    te <- try(fisher.test(tab), silent = TRUE)
    if (!inherits(te, "try-error")) {
      enrich_list[[i]] <- data.frame(
        var1 = pb_var$pb[i],
        var2 = pb_var$test[i],
        p_value = te$p.value,
        odds_ratio = as.numeric(te$estimate),
        total_expressed = nrow(sc_subset),
        degs = sum(sc_subset$sig),
        kegg_true = sum(sc_subset$sen_gene & sc_subset$sig),
        total_expressed_kegg = sum(sc_subset$sen_gene)
      )
    }
  }
}
enrich_df <- do.call(rbind, enrich_list)
print(enrich_df[order(match(enrich_df$var2, percentages)), ])
# saveRDS(enrich_df, enrich_rds)

# -------- Fisher test: AUCell labels vs. ground truth --------
AUCX2 <- readRDS(aux_rds)

thresholds <- data.frame(
  q = c(0.70, 0.80, 0.90, 0.95, 0.99),
  label = c("30percent", "20percent", "10percent", "5percent", "1percent"),
  stringsAsFactors = FALSE
)

meta_df <- readRDS(meta_rds)    # has exp, barcode, Sample_ID 

fisher_rows <- vector("list", nrow(thresholds))
for (x in seq_len(nrow(thresholds))) {
  AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
  AUC_mat$cell <- rownames(AUC_mat)
  thr <- quantile(AUC_mat$sen_list, thresholds$q[x])
  AUC_mat$sen <- AUC_mat$sen_list >= thr

  allm <- merge(AUC_mat, meta_df, by.x = "cell", by.y = "barcode")
  allm$category <- ifelse(allm$exp == "WT", "senescent", "control")

  te <- fisher.test(table(allm$sen, allm$category), alternative = "greater")
  fisher_rows[[x]] <- data.frame(
    test = thresholds$label[x],
    odds_ratio = as.numeric(te$estimate),
    p_value = te$p.value
  )
}
fisher_df <- do.call(rbind, fisher_rows)
print(fisher_df[order(match(fisher_df$test, percentages)), ])

# Optional quick correlation summary (AUCell vs truth)
AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
AUC_mat$cell <- rownames(AUC_mat)
allm <- merge(AUC_mat, meta_df, by.x = "cell", by.y = "barcode")
allm$sen_truth <- as.integer(allm$exp == "WT")
cat("Spearman AUCell vs WT truth:", cor(allm$sen_list, allm$sen_truth, method = "spearman"), "\n")

# -------- TRUE senescence PSEUDOBULK & DE --------
# build boolean and Sample_ID (reuse meta_df)
meta_true <- meta_df
meta_true$sen <- meta_true$exp == "WT"
if (!"Sample_ID" %in% names(meta_true) || any(is.na(meta_true$Sample_ID))) {
  meta_true$Sample_ID <- str_extract(meta_true$barcode, "(?<=-)[^-_]+$")
  meta_true$Sample_ID[is.na(meta_true$Sample_ID) | meta_true$Sample_ID == ""] <- as.character(meta_true$exp)
}

# attach to sce & aggregate
colData(sce)$sen <- meta_true$sen[match(colnames(sce), meta_true$barcode)]
colData(sce)$Sample_ID <- meta_true$Sample_ID[match(colnames(sce), meta_true$barcode)]

sene_df <- data.frame(
  cell     = colnames(sce),
  sen      = colData(sce)$sen,
  Sample_ID= colData(sce)$Sample_ID,
  stringsAsFactors = FALSE
)
sene_df$type <- "lung"
sene_df$cell_sen <- paste0(sene_df$type, "_", sene_df$sen)

sce_true <- sce
colData(sce_true)$sen_auc   <- sene_df$cell_sen[match(colnames(sce_true), sene_df$cell)]
colData(sce_true)$Sample_ID <- sene_df$Sample_ID[match(colnames(sce_true), sene_df$cell)]

pb_true <- aggregateToPseudoBulk(
  sce_true, assay = assay_name,
  cluster_id = "sen_auc",
  sample_id  = "Sample_ID",
  BPPARAM = BiocParallel::MulticoreParam(30)
)
# saveRDS(pb_true, file.path(true_expr_dir, "pseudobulk_truesen.RDS"))

# DE for true senescence
pb <- readRDS(file.path(true_expr_dir, "pseudobulk_truesen.RDS"))

# add ncells for dreamlet metadata
ncells  <- as.data.frame(int_colData(pb)$n_cells)
ncells$celltype <- rownames(ncells)
ncells2 <- ncells %>% pivot_longer(!celltype, names_to = "id", values_to = "count") %>% as.data.frame()
ncells2$id <- gsub("\\.", "_", ncells2$id)
colnames(ncells2) <- c("sen_auc", "Sample_ID", "ncell")

metadata(pb)$aggr_means$Indv_ID <- metadata(pb)$aggr_means$Sample_ID
metadata(pb)$aggr_means$Indv_ID <- gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
metadata(pb)$aggr_means$ncells  <- as.numeric(ncells2$ncell)

ct.pairs <- c("lung_TRUE","lung_FALSE")
fit_true <- try(dreamletCompareClusters_edgeR(
  pb, ct.pairs,
  min.samples = 1, min.cells = 1, min.count = 1,
  method = "none", min.prop = 0.4,
  useProcessOneAssay_edgeR = TRUE
), silent = TRUE)

if (!inherits(fit_true, "try-error")) {
  fit_top_true <- as.data.frame(topTags(fit_true, n = Inf))
  # saveRDS(fit_top_true, file.path(true_de_dir, "correct_edger_true_sen_fit.RDS"))
}

# -------- Compare predicted DE vs true-sen DE --------
true_sen <- readRDS(file.path(true_de_dir, "correct_edger_true_sen_fit.RDS"))
true_sen$symbol <- rownames(true_sen)

pred_list <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_dir, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(f)) next
  df <- readRDS(f)
  df$test   <- pb_var$test[i]
  df$symbol <- rownames(df)
  pred_list[[i]] <- df
}
predicted_sen <- do.call(rbind, pred_list)

merged <- merge(true_sen, predicted_sen, by = "symbol", all.y = TRUE)
targets <- c("1percent","5percent","10percent","20percent","30percent")
cor_result <- lapply(targets, function(lbl) {
  df <- subset(merged, test == lbl)
  if (nrow(df) == 0) return(data.frame(test = lbl, rho = NA, pval = NA))
  z <- suppressWarnings(cor.test(df$logFC.x, df$logFC.y, method = "spearman"))
  data.frame(test = lbl, rho = as.numeric(z$estimate), pval = z$p.value)
})
cor_result <- do.call(rbind, cor_result)
print(cor_result)
