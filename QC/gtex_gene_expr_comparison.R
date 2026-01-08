run_gtex_brain_corr <- function(
  # ---- Inputs ----
  gtex_counts_gct_path,
  gtex_sample_attr_ds_path,
  gtex_subject_pheno_ds_path,
  all_degs_rds_path = NULL,        
  degs_filter_rows  = c("exc_5percent","mg_30percent","int_5percent","ast_1percent","oli_10percent","opc_30percent"),
  
  # ---- Outputs ----
  output_dir,
  cache_dir,
  hist_all_pdf   = file.path(output_dir, "hist_total_signif_pairs.pdf"),
  hist_dir_pdf   = file.path(output_dir, "hist_total_signif_pairs_by_dir.pdf"),
  heatmap_pdf    = file.path(output_dir, "gtex_sig_pos_degs.pdf"),
  res_final_rds  = file.path(cache_dir,  "resFinal_brain_corr.RDS"),
  pairwise_tsv   = file.path(output_dir, "pairwise_breakdown.tsv"),
  gene_counts_tsv= file.path(output_dir, "gene_signif_counts.tsv"),
  
  # ---- Parameters ----
  cpm_thresh = 1,
  expr_prop  = 0.10,  
  fdr_alpha  = 0.05,
  min_common_subjects = 3,   
  cores = 40,
  seed  = 1,
  force_recompute = FALSE     
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(limma)
    library(assertthat)
    library(foreach)
    library(doParallel)
    library(ggplot2)
  })
  
  # ---- setup ----
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir,  showWarnings = FALSE, recursive = TRUE)
  
  must <- function(p) if (!file.exists(p)) stop("Missing file: ", p)
  lapply(list(gtex_counts_gct_path, gtex_sample_attr_ds_path, gtex_subject_pheno_ds_path), must)
  if (!is.null(all_degs_rds_path)) must(all_degs_rds_path)
  
  set.seed(seed)
  
  # ---------- Load GTEx counts ----------
  geneexp_all <- fread(gtex_counts_gct_path, data.table = FALSE)
  stopifnot(all(c("Name","Description") %in% colnames(geneexp_all)))
  rownames(geneexp_all) <- geneexp_all$Name
  geneexp_all$Name <- NULL
  geneexp_all$Description <- NULL
  
  # ---------- Metadata ----------
  sampleMeta <- fread(gtex_sample_attr_ds_path, data.table = FALSE)
  subjectMeta <- fread(gtex_subject_pheno_ds_path, data.table = FALSE)
  sampleMeta$SUBJID <- vapply(strsplit(sampleMeta$SAMPID, "-", fixed = TRUE),
                              function(x) paste0(x[1], "-", x[2]), character(1))
  sampleMeta <- merge(sampleMeta, subjectMeta, by = "SUBJID")
  keep <- intersect(sampleMeta$SAMPID, colnames(geneexp_all))
  geneexp_all <- geneexp_all[, keep, drop = FALSE]
  sampleMeta  <- sampleMeta[match(colnames(geneexp_all), sampleMeta$SAMPID), , drop = FALSE]
  assert_that(identical(colnames(geneexp_all), sampleMeta$SAMPID))
  
  # ---------- Brain subset ----------
  brainMeta <- subset(sampleMeta, SMTS == "Brain")
  brainExp  <- geneexp_all[, brainMeta$SAMPID, drop = FALSE]
  assert_that(identical(colnames(brainExp), brainMeta$SAMPID))
  
  # ---------- Expressed genes + voom ----------
  isexpr <- rowSums(cpm(brainExp) >= cpm_thresh) >= expr_prop * ncol(brainExp)
  genesAll <- DGEList(counts = brainExp[isexpr, ])
  genesAll <- calcNormFactors(genesAll)
  vobj <- limma::voom(genesAll, plot = FALSE)
  
  # ---------- All brain region pairs ----------
  regions <- unique(brainMeta$SMTSD)
  combs <- t(combn(regions, 2))
  
  # ---------- Cache heavy correlation results ----------
  need_compute <- force_recompute || !file.exists(res_final_rds)
  if (need_compute) {
    registerDoParallel(cores = cores)
    final <- foreach(comb = 1:nrow(combs), .combine = rbind) %dopar% {
      region1 <- combs[comb, 1]; region2 <- combs[comb, 2]
      meta1 <- brainMeta[brainMeta$SMTSD == region1, ]
      meta2 <- brainMeta[brainMeta$SMTSD == region2, ]
      commonSubjects <- intersect(meta1$SUBJID, meta2$SUBJID)
      # guard: skip pairs with too few overlapping subjects
      if (length(commonSubjects) < min_common_subjects) {
        return(data.frame(region1 = character(0), region2 = character(0),
                          gene = character(0), spearmanCor = numeric(0), pval = numeric(0)))
      }
      meta1 <- meta1[match(commonSubjects, meta1$SUBJID), ]
      meta2 <- meta2[match(commonSubjects, meta2$SUBJID), ]
      assert_that(identical(meta1$SUBJID, meta2$SUBJID))
      v1 <- vobj[, meta1$SAMPID]; v2 <- vobj[, meta2$SAMPID]
      assert_that(identical(colnames(v1), meta1$SAMPID))
      assert_that(identical(colnames(v2), meta2$SAMPID))
      
      genes <- intersect(rownames(v1), rownames(v2))
      res_list <- lapply(genes, function(g) {
        ct <- suppressWarnings(cor.test(v1$E[g, ], v2$E[g, ], method = "spearman"))
        data.frame(region1 = region1, region2 = region2, gene = g,
                   spearmanCor = unname(ct$estimate), pval = ct$p.value,
                   stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, res_list)
      as.data.frame(out, stringsAsFactors = FALSE)
    }
    resFinal <- as.data.frame(final, stringsAsFactors = FALSE)
    if (nrow(resFinal) > 0) {
      resFinal$adj.p.val <- p.adjust(resFinal$pval, method = "fdr")
    } else {
      resFinal$adj.p.val <- numeric(0)
    }
    saveRDS(resFinal, res_final_rds)
  } else {
    resFinal <- readRDS(res_final_rds)
    resFinal <- as.data.frame(resFinal, stringsAsFactors = FALSE)
  }
  
  # ---------- Gene-level counts (how many significant pairs) ----------
  registerDoParallel(cores = cores)
  unique_genes <- unique(resFinal$gene)
  gene_counts <- foreach(g = unique_genes, .combine = rbind) %dopar% {
    sub <- resFinal[resFinal$gene == g, ]
    data.frame(
      gene = g,
      total_signif_pairs = sum(sub$adj.p.val < fdr_alpha),
      total_signif_pairs_positive = sum(sub$adj.p.val < fdr_alpha & sub$spearmanCor > 0),
      total_signif_pairs_negative = sum(sub$adj.p.val < fdr_alpha & sub$spearmanCor < 0),
      stringsAsFactors = FALSE
    )
  }
  gene_counts <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
  fwrite(gene_counts, gene_counts_tsv, sep = "\t")
  
  # ---------- Pairwise breakdown ----------
  pairwise <- foreach(i = 1:nrow(combs), .combine = rbind) %dopar% {
    r1 <- combs[i, 1]; r2 <- combs[i, 2]
    sub <- resFinal[resFinal$region1 == r1 & resFinal$region2 == r2, ]
    nS  <- sum(sub$adj.p.val < fdr_alpha)
    data.frame(
      region1 = r1, region2 = r2,
      numSignif = nS,
      numSignifPos = sum(sub$adj.p.val < fdr_alpha & sub$spearmanCor > 0),
      numSignifNeg = sum(sub$adj.p.val < fdr_alpha & sub$spearmanCor < 0),
      numNotSignif = sum(sub$adj.p.val >= fdr_alpha),
      stringsAsFactors = FALSE
    )
  }
  pairwise <- as.data.frame(pairwise, stringsAsFactors = FALSE)
  pairwise$proportionSignif    <- with(pairwise, ifelse(numSignif + numNotSignif > 0, numSignif / (numSignif + numNotSignif), NA_real_))
  pairwise$proportionSignifPos <- with(pairwise, ifelse(numSignif + numNotSignif > 0, numSignifPos / (numNotSignif + numSignif), NA_real_))
  
  # ---------- Optional DE overlay (Anina lists) ----------
  if (!is.null(all_degs_rds_path)) {
    all_degs_test <- readRDS(all_degs_rds_path)
    sub_de <- all_degs_test[ all_degs_test$de_test %in% c("sen","bulk_no_odc","single"), ]
    sub_de <- sub_de[ sub_de$de_test != "sen" | (sub_de$de_test == "sen" & sub_de$feature %in% degs_filter_rows), ]
    resFinal$ensemblNoVersion <- sub("\\..*$", "", resFinal$gene)
    de_ids <- unique(sub_de$ID[sub_de$ID != ""])
    resFinal$is_DE <- resFinal$ensemblNoVersion %in% de_ids
    total_DE <- length(unique(resFinal$ensemblNoVersion[resFinal$is_DE]))
    pairwise$numSignifPosAndDE <- vapply(1:nrow(pairwise), function(i) {
      r1 <- pairwise$region1[i]; r2 <- pairwise$region2[i]
      s <- resFinal[resFinal$region1 == r1 & resFinal$region2 == r2, ]
      sum(s$adj.p.val <= fdr_alpha & s$spearmanCor > 0 & s$is_DE)
    }, integer(1))
    pairwise$proportionSignifPosAndDE <- if (total_DE > 0) pairwise$numSignifPosAndDE / total_DE else NA_real_
  }
  
  fwrite(pairwise, pairwise_tsv, sep = "\t")
  
  # ---------- Plots ----------
  # histograms
  gcounts_posneg <- {
    resCountPos <- gene_counts[, c("gene","total_signif_pairs_positive")]
    names(resCountPos) <- c("gene","total_signif_pairs"); resCountPos$direction <- "positive"
    resCountNeg <- gene_counts[, c("gene","total_signif_pairs_negative")]
    names(resCountNeg) <- c("gene","total_signif_pairs"); resCountNeg$direction <- "negative"
    rbind(resCountPos, resCountNeg)
  }
  p1 <- ggplot(gene_counts, aes(total_signif_pairs)) + geom_histogram(bins = 78) +
        theme_bw() + labs(x = "Significant region-pairs per gene", y = "Count")
  ggsave(hist_all_pdf, p1, width = 7, height = 5)
  
  p2 <- ggplot(gcounts_posneg, aes(total_signif_pairs)) + 
        geom_histogram(bins = 78) + facet_wrap(~direction, scales = "free") +
        theme_bw() + labs(x = "Significant region-pairs per gene", y = "Count")
  ggsave(hist_dir_pdf, p2, width = 9, height = 4.5)
  
  # heatmap (if DE overlay available, use proportionSignifPosAndDE; else proportionSignifPos)
  heat_fill <- if (!is.null(all_degs_rds_path)) "proportionSignifPosAndDE" else "proportionSignifPos"
  
  # Pretty labels like your script
  pretty_lab <- function(x) {
    x <- sub("^Brain - ", "", x)
    x <- gsub("·", " ", x)
    x <- gsub("Nucleus accumbens \\(basal ganglia\\)", "Nucleus accumbens\n(basal ganglia)", x)
    x <- gsub("Anterior cingulate cortex \\(BA24\\)", "Anterior cingulate\ncortex (BA24)", x)
    x
  }
  pairwise$region1_lab <- pretty_lab(pairwise$region1)
  pairwise$region2_lab <- pretty_lab(pairwise$region2)
  
  # fix ordering to match unique combs (prettified)
  combs_lab <- cbind(pretty_lab(combs[,1]), pretty_lab(combs[,2]))
  pairwise$region1_lab <- factor(pairwise$region1_lab, unique(combs_lab[,1]))
  pairwise$region2_lab <- factor(pairwise$region2_lab, unique(combs_lab[,2]))
  
  d <- ggplot(pairwise, aes(x = region1_lab, y = region2_lab, fill = .data[[heat_fill]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[heat_fill]])), size = 4, color = "black") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 11),
      axis.ticks  = element_line(color = "black"),
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.box.just = "right",
      legend.title.align = 0.5,
      legend.text.align = 0.5,
      legend.margin = margin(6, 6, 6, 6),
      legend.key = element_rect(color = "black", size = 0.5),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    labs(x = "Brain Region", y = "Brain Region",
         fill = if (!is.null(all_degs_rds_path))
           "Proportion of DE genes with\nSignificant Positive Correlations"
         else
           "Proportion of genes with\nSignificant Positive Correlations") +
    guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.position = "top"))
  
  ggsave(heatmap_pdf, d, width = 10, height = 9)
  
  invisible(list(
    n_regions = length(regions),
    n_pairs = nrow(combs),
    n_genes = nrow(vobj$E),
    hist_all_pdf = hist_all_pdf,
    hist_dir_pdf = hist_dir_pdf,
    heatmap_pdf = heatmap_pdf,
    res_final_rds = res_final_rds,
    pairwise_tsv = pairwise_tsv,
    gene_counts_tsv = gene_counts_tsv
  ))
}


res <- run_gtex_brain_corr(
  # inputs
  gtex_counts_gct_path       = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
  gtex_sample_attr_ds_path   = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
  gtex_subject_pheno_ds_path = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
  all_degs_rds_path          = "all_degs_test.RDS", 
  
  # outputs
  output_dir = "/plots/",                     
  cache_dir  = "/tmp",                
  heatmap_pdf = "gtex_sig_pos_degs.pdf", 
  
  # params
  cpm_thresh = 1,
  expr_prop  = 0.10,
  fdr_alpha  = 0.05,
  cores = 40,
  seed  = 1
)
