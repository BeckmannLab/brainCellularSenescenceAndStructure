# -------- DATA DOWNLOAD (shell) --------
# wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115301/suppl/GSE115301%5FGrowing%5FSen%5F10x%5Fcount%2Etxt%2Egz"   # counts
# wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115301/suppl/GSE115301%5FGrowing%5FSen%5F10x%5Fmetadata%2Etxt%2Egz" # metadata

# -------- LIBRARIES --------
options(stringsAsFactors = FALSE)

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
library(tidyr)
library(stringr)
library(doParallel)
registerDoParallel(cores = 10)

# Custom helper functions (edgeR wrappers)
source("dreamletCompareClusters_edgeR.R")
source("processOneAssay_edgeR.R")
source("dreamletCompareClusters_edgeR_notPaired.R")

# -------- FILE PATHS --------
base            <- "GSE115301/"
aucell_out      <- file.path(base, "aucell/")
path_results    <- file.path(aucell_out, "cutoffs/")
expr_results    <- file.path(aucell_out, "expression/")
de_out          <- file.path(aucell_out, "DE/")
true_sen_folder <- file.path(aucell_out, "true_senescenceDE/")
true_sen_expression <- file.path(true_sen_folder, "expression/")
true_sen_results    <- file.path(true_sen_folder, "DE/")

# outputs
enrichment_results <- file.path(aucell_out, "enrichment_or_12.16.25.RDS")
meta_file_save     <- file.path(base, "meta_for_counts.RDS")
sce_save           <- file.path(base, "sc_count_with_meta.RDS")
aux_save           <- file.path(base, "AUCell_GSE115301_1.16.24.RDS")
auc_mat_path       <- file.path(path_results, "all_lung_AUC_mat.RDS")

# ensure directories exist
invisible(lapply(
  c(aucell_out, path_results, expr_results, de_out,
    true_sen_folder, true_sen_expression, true_sen_results),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# -------- LOAD COUNTS & METADATA --------
meta  <- fread(file.path(base, "GSE115301_Growing_Sen_10x_metadata.txt"))
count <- fread(file.path(base, "GSE115301_Growing_Sen_10x_count.txt"))

# counts: first column is gene symbol, remaining columns are cell barcodes
count <- as.data.frame(count)
rownames(count) <- count[[1]]
count2 <- count[, -1, drop = FALSE]

# metadata: first col = barcode (V1), second col = condition (V2)
meta2 <- as.data.frame(meta)
colnames(meta2)[1:2] <- c("barcode", "group")
meta2 <- meta2[!duplicated(meta2$barcode), ]
rownames(meta2) <- meta2$barcode

# -------- ALIGN META TO COUNTS (fixes SCE dim error) --------
cells_counts <- colnames(count2)
keep <- intersect(cells_counts, rownames(meta2))
stopifnot(length(keep) > 0)

count2 <- count2[, keep, drop = FALSE]
meta2  <- meta2 [ keep,  , drop = FALSE]

# convenience columns (status/type/Sample_ID)
meta2$status    <- ifelse(meta2$group == "Senescence1", "senescence", "control")
meta2$type      <- "lung"
meta2$Sample_ID <- str_extract(rownames(meta2), "(?<=-).*")

saveRDS(meta2, meta_file_save)

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
all_genes <- unique(c(sen_genes, lung_genes))
stopifnot(length(intersect(rownames(count), all_genes)) > 0)
geneSets <- GeneSet(all_genes, setName = "sen_list")

# -------- BUILD SCE --------
countMat <- as(as.matrix(count2), "dgCMatrix")
sce <- SingleCellExperiment(assays = list(counts = countMat), colData = meta2)
saveRDS(sce, sce_save)

# -------- RUN AUCELL --------
AUCX2 <- AUCell_run(countMat, geneSets = geneSets, BPPARAM = BiocParallel::MulticoreParam(20))
saveRDS(AUCX2, aux_save)

# -------- THRESHOLD CELLS BY AUCELL SCORE AT MULTIPLE QUANTILES --------
percent_labels <- c("1percent","5percent","10percent","20percent","30percent")
quantiles      <- c(0.99, 0.95, 0.90, 0.80, 0.70)
all_var <- data.frame(percent = percent_labels, q = quantiles)

final_list <- vector("list", nrow(all_var))
for (i in seq_len(nrow(all_var))) {
  AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
  AUC_mat$cell <- rownames(AUC_mat)
  thr <- quantile(AUC_mat$sen_list, all_var$q[i])
  AUC_mat$sen <- AUC_mat$sen_list >= thr
  AUC_mat$percent <- all_var$percent[i]
  q_lbl <- gsub("\\.", "_", sprintf("%.2f", all_var$q[i]))
  saveRDS(AUC_mat, file.path(path_results, sprintf("%s_lung_AUC_mat_top_%s.RDS", q_lbl, all_var$percent[i])))
  final_list[[i]] <- AUC_mat
}
all_sene <- do.call(rbind, final_list)
saveRDS(all_sene, auc_mat_path)

# -------- PSEUDOBULK (predicted labels) --------
all_sene <- readRDS(auc_mat_path)
all_sene$type <- "lung"

percentages <- percent_labels
pb_var <- data.frame(
  file  = paste0("pseudobulk_", percentages, ".RDS"),
  test  = percentages,
  pb    = paste0("pb_", percentages),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pb_var))) {
  sene <- subset(all_sene, percent == pb_var$test[i])
  sene$cell_sen <- paste0(sene$type, "_", sene$sen)
  sen3 <- sene[, c("cell_sen","cell"), drop = FALSE]

  matches <- match(colnames(countMat), sen3$cell)
  valid   <- !is.na(matches)
  sce2    <- sce[, valid]
  colData(sce2)$sen_auc   <- as.character(sen3[matches[valid], "cell_sen", drop = TRUE])
  colData(sce2)$Sample_ID <- colData(sce)[, "Sample_ID"][valid]

  out_file <- file.path(expr_results, pb_var$file[i])
  if (!file.exists(out_file)) {
    pb <- aggregateToPseudoBulk(
      sce2, assay = "counts",
      cluster_id = "sen_auc",
      sample_id  = "Sample_ID",
      BPPARAM = BiocParallel::MulticoreParam(30)
    )
    saveRDS(pb, out_file)
  }
}

# -------- DIFFERENTIAL EXPRESSION (edgeR via dreamlet wrapper) --------
for (i in seq_len(nrow(pb_var))) {
  in_file  <- file.path(expr_results, pb_var$file[i])
  out_file <- file.path(de_out, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(in_file) || file.exists(out_file)) next

  pb <- readRDS(in_file)

  # add ncells to metadata
  ncells  <- as.data.frame(int_colData(pb)$n_cells)
  ncells$celltype <- rownames(ncells)
  ncells2 <- ncells %>%
    pivot_longer(!celltype, names_to = "id", values_to = "count") %>%
    as.data.frame()
  ncells2$id <- gsub("\\.", "_", ncells2$id)
  colnames(ncells2) <- c("sen_auc", "Sample_ID", "ncell")

  metadata(pb)$aggr_means$Indv_ID <- metadata(pb)$aggr_means$Sample_ID
  metadata(pb)$aggr_means$Indv_ID <- gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
  metadata(pb)$aggr_means$ncells  <- as.numeric(ncells2$ncell)

  ct.pairs <- c("lung_TRUE", "lung_FALSE")
  fit <- try(dreamletCompareClusters_edgeR(
    pb, ct.pairs,
    min.samples = 1, min.cells = 1, min.count = 1,
    method = "fixed", min.prop = 0.4,
    useProcessOneAssay_edgeR = TRUE
  ), silent = TRUE)

  if (!inherits(fit, "try-error")) {
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    saveRDS(fit_top, out_file)
  }
}

# -------- SUMMARIZE DE RESULTS --------
summary_list <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_out, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(f)) next
  df <- readRDS(f)
  summary_list[[i]] <- data.frame(
    sig              = sum(df$FDR <= 0.05),
    total            = nrow(df),
    percentage       = sum(df$FDR <= 0.05) / nrow(df),
    celltype         = "lung",
    test             = pb_var$test[i],
    percentage_up    = sum(df$logFC > 0) / nrow(df),
    percentage_down  = sum(df$logFC < 0) / nrow(df),
    pi1              = 1 - pi0est(df$PValue)$pi0
  )
}
summary_df <- do.call(rbind, summary_list)
print(summary_df[order(match(summary_df$test, percentages)), ])

# -------- Enrichment of senescence genes among DE --------
all_de <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_out, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
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

results <- list()
for (i in seq_len(nrow(pb_var))) {
  sc_subset <- subset(sc_with_sen, test == pb_var$test[i])
  tab <- table(sc_subset$sen_gene, sc_subset$sig)
  if (nrow(tab) >= 2 && ncol(tab) >= 2) {
    te <- try(fisher.test(tab), silent = TRUE)
    if (!inherits(te, "try-error")) {
      results[[i]] <- data.frame(
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
enrich_df <- do.call(rbind, results)
print(enrich_df[order(match(enrich_df$var2, percentages)), ])
saveRDS(enrich_df, enrichment_results)

# -------- Fisher test: AUCell labels vs. ground truth --------
AUCX2 <- readRDS(aux_save)
thresholds <- data.frame(
  q = c(0.70, 0.80, 0.90, 0.95, 0.99),
  label = c("30percent", "20percent", "10percent", "5percent", "1percent"),
  stringsAsFactors = FALSE
)

fisher_list <- list()
for (x in seq_len(nrow(thresholds))) {
  AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
  AUC_mat$cell <- rownames(AUC_mat)
  thr <- quantile(AUC_mat$sen_list, thresholds$q[x])
  AUC_mat$sen <- AUC_mat$sen_list >= thr

  allm <- merge(AUC_mat, meta2, by.x = "cell", by.y = "barcode")
  allm$category <- ifelse(allm$group == "Senescence1", "sen", "control")

  te <- fisher.test(table(allm$sen, allm$category), alternative = "greater")
  fisher_list[[x]] <- data.frame(
    test = thresholds$label[x],
    odds_ratio = as.numeric(te$estimate),
    p_value = te$p.value
  )
}
fisher_df <- do.call(rbind, fisher_list)
print(fisher_df[order(match(fisher_df$test, percentages)), ])

# -------- Correlate AUCell scores with truth --------
AUC_mat <- as.data.frame(t(getAUC(AUCX2)))
AUC_mat$cell <- rownames(AUC_mat)
allm <- merge(AUC_mat, meta2, by.x = "cell", by.y = "barcode")
allm$sen_truth <- as.integer(allm$group == "Senescence1")
cat("Spearman AUCell vs truth:", cor(allm$sen_list, allm$sen_truth, method = "spearman"), "\n")

# -------- TRUE senescence PSEUDOBULK & DE (robust condition detection) --------
sce <- readRDS(sce_save)
countMat <- assay(sce)

meta_true <- as.data.frame(colData(sce))

# pick first condition-like column that exists
cond_candidates <- c("group", "V2", "status", "condition", "Condition", "Group")
cond_col <- intersect(cond_candidates, colnames(meta_true))
if (length(cond_col) == 0) {
  stop("No condition column found in colData(sce). Have: ",
       paste(colnames(meta_true), collapse = ", "))
}
cond_col <- cond_col[1]

vals <- tolower(as.character(meta_true[[cond_col]]))
meta_true$sen <- vals %in% c("senescence1","senescence","sen","true","1","yes")

# ensure Sample_ID in colData
if (!"Sample_ID" %in% names(meta_true) || any(is.na(meta_true$Sample_ID))) {
  meta_true$Sample_ID <- if ("Sample_ID" %in% colnames(colData(sce))) {
    as.character(colData(sce)$Sample_ID)
  } else {
    str_extract(rownames(meta_true), "(?<=-).*")
  }
}
colData(sce)$sen <- meta_true$sen
colData(sce)$Sample_ID <- meta_true$Sample_ID

sene <- meta_true
sene$cell <- rownames(sene)
sene$type <- "lung"
sene$cell_sen <- paste0(sene$type, "_", sene$sen)
sen3 <- sene[, c("cell_sen","cell","Sample_ID"), drop = FALSE]

matches <- match(colnames(countMat), sen3$cell)
valid   <- !is.na(matches)
sce2    <- sce[, valid]
colData(sce2)$sen_auc   <- as.character(sen3[matches[valid], "cell_sen", drop = TRUE])
colData(sce2)$Sample_ID <- as.character(sen3[matches[valid], "Sample_ID", drop = TRUE])

pb_true <- aggregateToPseudoBulk(
  sce2, assay = "counts",
  cluster_id = "sen_auc",
  sample_id  = "Sample_ID",
  BPPARAM = BiocParallel::MulticoreParam(30)
)
# saveRDS(pb_true, file.path(true_sen_expression, "pseudobulk_truesen.RDS"))

# DE for true senescence
pb <- readRDS(file.path(true_sen_expression, "pseudobulk_truesen.RDS"))

# add ncells
ncells  <- as.data.frame(int_colData(pb)$n_cells)
ncells$celltype <- rownames(ncells)
ncells2 <- ncells %>% pivot_longer(!celltype, names_to = "id", values_to = "count") %>% as.data.frame()
ncells2$id <- gsub("\\.", "_", ncells2$id)
colnames(ncells2) <- c("sen_auc", "Sample_ID", "ncell")

metadata(pb)$aggr_means$Indv_ID <- metadata(pb)$aggr_means$Sample_ID
metadata(pb)$aggr_means$Indv_ID <- gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
metadata(pb)$aggr_means$ncells  <- as.numeric(ncells2$ncell)

ct.pairs <- c("lung_TRUE", "lung_FALSE")
fit_true <- try(dreamletCompareClusters_edgeR(
  pb, ct.pairs,
  min.samples = 1, min.cells = 1, min.count = 1,
  method = "none", min.prop = 0.4,
  useProcessOneAssay_edgeR = TRUE
), silent = TRUE)

if (!inherits(fit_true, "try-error")) {
  fit_top_true <- as.data.frame(topTags(fit_true, n = Inf))
  # saveRDS(fit_top_true, file.path(true_sen_results, "correct_edger_true_sen_fit.RDS"))
}

# -------- Compare predicted DE vs true-sen DE (Spearman on logFC) --------
true_sen <- readRDS(file.path(true_sen_results, "correct_edger_true_sen_fit.RDS"))
true_sen$symbol <- rownames(true_sen)

pred_list <- list()
for (i in seq_len(nrow(pb_var))) {
  f <- file.path(de_out, sprintf("correct_edger_%s_fit.RDS", pb_var$test[i]))
  if (!file.exists(f)) next
  df <- readRDS(f)
  df$test   <- pb_var$test[i]
  df$symbol <- rownames(df)
  pred_list[[i]] <- df
}
predicted_sen <- do.call(rbind, pred_list)

merged <- merge(true_sen, predicted_sen, by = "symbol", all.y = TRUE)
targets <- c("10percent", "20percent", "30percent")
cor_result <- lapply(targets, function(lbl) {
  df <- subset(merged, test == lbl)
  if (nrow(df) == 0) return(data.frame(test = lbl, rho = NA, pval = NA))
  z <- suppressWarnings(cor.test(df$logFC.x, df$logFC.y, method = "spearman"))
  data.frame(test = lbl, rho = as.numeric(z$estimate), pval = z$p.value)
})
cor_result <- do.call(rbind, cor_result)
print(cor_result)
