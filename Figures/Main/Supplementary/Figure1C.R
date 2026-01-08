

run_gtex_pfc_mds <- function(
  # -------- Input files --------
  rosmap_counts_path,
  msbb_counts_path,
  lbp_counts_path,
  gtex_counts_gct_path,
  gtex_sample_attr_ds_path,
  gtex_subject_pheno_ds_path,
  # Dictionaries are optional (unused in downstream steps, kept for parity)
  gtex_sample_attr_dd_xlsx = NULL,
  gtex_subject_pheno_dd_xlsx = NULL,
  
  # -------- Output / cache --------
  output_dir,
  cache_dir,
  xlsx_out = file.path(output_dir, "mds_pfc_other_tissues.xlsx"),
  pdf_out_alt  = file.path(output_dir, "TMP_plot2.pdf"),
  
  # -------- Parameters --------
  read_lengths = list(rosmap = 101, msbb = 100, lbp = 100, gtex = 76),
  tpm_threshold = 0.5,         # gene is "expressed" if TPM>= threshold in >= expr_prop * samples
  expr_prop = 0.10,
  threads = 30,
  seed = 123
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(assertthat)
    library(parallelDist)
    library(ggplot2)
    library(writexl)
  })
  
  # ---- helpers --------------------------------------------------------------
  must_exist <- function(p) if (!file.exists(p)) stop("File not found: ", p)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir,  showWarnings = FALSE, recursive = TRUE)
  
  # robust TPM (same math you used; vectorized & with checks)
  counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
    stopifnot(is.matrix(counts) || is.data.frame(counts))
    counts <- as.matrix(counts)
    stopifnot(length(featureLength) == nrow(counts))
    stopifnot(length(meanFragmentLength) == ncol(counts))
    
    effLen <- sweep(matrix(featureLength, nrow = length(featureLength), ncol = ncol(counts)),
                    2, meanFragmentLength, FUN = function(fl, mfl) fl - mfl + 1)
    keep <- apply(effLen, 1, function(x) min(x) > 1)
    counts <- counts[keep, , drop = FALSE]
    effLen <- effLen[keep, , drop = FALSE]
    featureLength <- featureLength[keep]
    
    # log-sum-exp per column
    rate <- log(pmax(counts, 0) + 1e-8) - log(effLen)  # tiny epsilon to avoid -Inf
    lse <- apply(rate, 2, function(col) {
      m <- max(col)
      m + log(sum(exp(col - m)))
    })
    tpm <- exp(sweep(rate, 2, lse, "-") + log(1e6))
    rownames(tpm) <- rownames(counts)
    colnames(tpm) <- colnames(counts)
    tpm
  }
  
  # ---- read inputs ----------------------------------------------------------
  lapply(c(rosmap_counts_path, msbb_counts_path, lbp_counts_path,
           gtex_counts_gct_path, gtex_sample_attr_ds_path, gtex_subject_pheno_ds_path), must_exist)
  
  message("Reading count matrices…")
  rosmap_raw <- fread(rosmap_counts_path, data.table = FALSE)
  msbb_raw   <- fread(msbb_counts_path,   data.table = FALSE)
  lbp_raw    <- fread(lbp_counts_path,    data.table = FALSE)
  
  # ROSMAP / MSBB assumed first 6 columns are bed-like, incl. Geneid (+ Length in MSBB)
  bed_rosmap <- rosmap_raw[, 1:6, drop = FALSE]
  bed_msbb   <- msbb_raw[,   1:6, drop = FALSE]
  if (!"Geneid" %in% colnames(bed_rosmap) || !"Geneid" %in% colnames(bed_msbb)) {
    stop("Expected a 'Geneid' column in the first 6 columns for ROSMAP/MSBB files.")
  }
  bed_rosmap$ensembl <- sub("\\..*$", "", bed_rosmap$Geneid)
  bed_msbb$ensembl   <- sub("\\..*$", "", bed_msbb$Geneid)
  
  # strip bed columns to counts
  rosmap_counts <- rosmap_raw[, 7:ncol(rosmap_raw), drop = FALSE]
  rownames(rosmap_counts) <- bed_rosmap$Geneid
  
  msbb_counts <- msbb_raw[, 7:ncol(msbb_raw), drop = FALSE]
  rownames(msbb_counts) <- bed_msbb$Geneid
  
  # LBP has Geneid in first col named V1
  lbp_counts <- lbp_raw
  if (!"V1" %in% colnames(lbp_counts)) stop("LBP counts expected a 'V1' column containing Geneid.")
  bed_lbp <- data.frame(Geneid = lbp_counts$V1, stringsAsFactors = FALSE)
  bed_lbp$ensembl <- sub("\\..*$", "", bed_lbp$Geneid)
  rownames(lbp_counts) <- bed_lbp$Geneid
  lbp_counts$V1 <- NULL
  
  message("Reading GTEx counts (GCT)…")
  gtex_raw <- fread(gtex_counts_gct_path, data.table = FALSE)
  if (!all(c("Name", "Description") %in% colnames(gtex_raw))) {
    stop("GTEx GCT expected columns 'Name' and 'Description'.")
  }
  rownames(gtex_raw) <- gtex_raw$Name
  gtex_raw$Name <- NULL
  gtex_raw$Description <- NULL
  bed_gtex <- data.frame(Geneid = rownames(gtex_raw), stringsAsFactors = FALSE)
  bed_gtex$ensembl <- sub("\\..*$", "", bed_gtex$Geneid)
  
  # ---- deduplicate by Ensembl (keep only genes mapping uniquely in each set) ----
  uniq <- function(x) x[x %in% names(which(table(x) == 1))]
  bed_msbb  <- bed_msbb [bed_msbb$ensembl  %in% uniq(bed_msbb$ensembl), ]
  bed_lbp   <- bed_lbp  [bed_lbp$ensembl   %in% uniq(bed_lbp$ensembl), ]
  bed_gtex  <- bed_gtex [bed_gtex$ensembl  %in% uniq(bed_gtex$ensembl), ]
  
  # Merge keys
  colnames(bed_gtex)[colnames(bed_gtex) == "Geneid"] <- "Geneid.gtex"
  bed_key <- merge(merge(bed_msbb[, c("Geneid", "ensembl")], 
                         bed_lbp[,  c("Geneid", "ensembl")], 
                         by = "ensembl",
                         suffixes = c(".msbb", ".lbp")),
                   bed_gtex, by = "ensembl")
  # For parity with your original naming
  colnames(bed_key)[colnames(bed_key) == "Geneid.msbb"] <- "Geneid.msbb.rosmap"
  
  # Subset/align counts by merged key
  rosmap_counts <- rosmap_counts[bed_key$Geneid.msbb.rosmap, , drop = FALSE]
  msbb_counts   <- msbb_counts  [bed_key$Geneid.msbb.rosmap, , drop = FALSE]
  lbp_counts    <- lbp_counts   [bed_key$Geneid.lbp,         , drop = FALSE]
  gtex_raw      <- gtex_raw     [bed_key$Geneid.gtex,        , drop = FALSE]
  
  rownames(rosmap_counts) <- bed_key$ensembl
  rownames(msbb_counts)   <- bed_key$ensembl
  rownames(lbp_counts)    <- bed_key$ensembl
  rownames(gtex_raw)      <- bed_key$ensembl
  
  # ---- Metadata (GTEx) ------------------------------------------------------
  message("Reading GTEx metadata…")
  sampleMeta <- fread(gtex_sample_attr_ds_path, data.table = FALSE)
  subjectMeta <- fread(gtex_subject_pheno_ds_path, data.table = FALSE)
  sampleMeta$SUBJID <- vapply(strsplit(sampleMeta$SAMPID, "-", fixed = TRUE),
                              function(x) paste0(x[1], "-", x[2]), character(1))
  sampleMeta <- merge(sampleMeta, subjectMeta, by = "SUBJID")
  
  # keep only GTEx columns we have
  gtex_samples <- intersect(sampleMeta$SAMPID, colnames(gtex_raw))
  gtex_raw <- gtex_raw[, gtex_samples, drop = FALSE]
  sampleMeta <- sampleMeta[match(colnames(gtex_raw), sampleMeta$SAMPID), , drop = FALSE]
  assert_that(identical(colnames(gtex_raw), sampleMeta$SAMPID))
  
  # ---- Build synthetic meta for other cohorts & assemble global counts -------
  message("Building combined metadata…")
  meta_cols <- c("SUBJID","SAMPID","SMATSSCR","SMCENTER","SMPTHNTS","SMRIN","SMTS","SMTSD",
                 "SMUBRID","SMTSISCH","SMTSPAX","SMNABTCH","SMNABTCHT","SMNABTCHD","SMGEBTCH",
                 "SMGEBTCHD","SMGEBTCHT","SMAFRZE","SMGTC","SME2MPRT","SMCHMPRS","SMNTRART",
                 "SMNUMGPS","SMMAPRT","SMEXNCRT","SM550NRM","SMGNSDTC","SMUNMPRT","SM350NRM",
                 "SMRDLGTH","SMMNCPB","SME1MMRT","SMSFLGTH","SMESTLBS","SMMPPD","SMNTERRT",
                 "SMRRNANM","SMRDTTL","SMVQCFL","SMMNCV","SMTRSCPT","SMMPPDPR","SMCGLGTH",
                 "SMGAPPCT","SMUNPDRD","SMNTRNRT","SMMPUNRT","SMEXPEFF","SMMPPDUN","SME2MMRT",
                 "SME2ANTI","SMALTALG","SME2SNSE","SMMFLGTH","SME1ANTI","SMSPLTRD","SMBSMMRT",
                 "SME1SNSE","SME1PCTS","SMRRNART","SME1MPRT","SMNUM5CD","SMDPMPRT","SME2PCTS",
                 "SEX","AGE","DTHHRDY")
  
  synth_meta <- function(counts, smts_label, read_len) {
    m <- as.data.frame(matrix(NA, nrow = ncol(counts), ncol = length(meta_cols)))
    colnames(m) <- meta_cols
    m$SUBJID <- colnames(counts)
    m$SAMPID <- colnames(counts)
    m$SMTS   <- smts_label
    m$SMTSD  <- smts_label
    m$readLength <- read_len
    m
  }
  
  meta_rosmap <- synth_meta(rosmap_counts, "PFC_rosmap", read_lengths$rosmap)
  meta_msbb   <- synth_meta(msbb_counts,   "PFC_msbb",   read_lengths$msbb)
  meta_lbp    <- synth_meta(lbp_counts,    "PFC_lbp",    read_lengths$lbp)
  sampleMeta$readLength <- read_lengths$gtex
  
  meta_all <- rbind(sampleMeta[, c(meta_cols, "readLength")],
                    meta_rosmap[, c(meta_cols, "readLength")],
                    meta_msbb[,   c(meta_cols, "readLength")],
                    meta_lbp[,    c(meta_cols, "readLength")])
  
  counts_all <- cbind(rosmap_counts, msbb_counts, lbp_counts, gtex_raw)
  # align meta
  meta_all <- meta_all[match(colnames(counts_all), meta_all$SAMPID), , drop = FALSE]
  assert_that(identical(colnames(counts_all), meta_all$SAMPID))
  assert_that(identical(rownames(counts_all), bed_key$ensembl))
  
  # ---- Feature lengths (from MSBB bed) --------------------------------------
  if (!"Length" %in% colnames(bed_msbb)) {
    stop("Expected a 'Length' column in MSBB bed (first 6 columns) to supply feature lengths.")
  }
  # Reorder to bed_key order
  feature_len <- bed_msbb$Length[match(bed_key$Geneid.msbb.rosmap, bed_msbb$Geneid)]
  names(feature_len) <- bed_key$ensembl
  
  # ---- TPM & QQ normalization (all genes) -----------------------------------
  message("Computing TPM (all genes)…")
  tpms_all <- counts_to_tpm(counts_all, 
                            featureLength = feature_len[rownames(counts_all)],
                            meanFragmentLength = meta_all$readLength)
  qq_all <- apply(tpms_all, 2, function(x) qqnorm(x, plot.it = FALSE)$x)
  rownames(qq_all) <- rownames(tpms_all)
  assert_that(identical(colnames(qq_all), meta_all$SAMPID))
  
  # ---- Distances / MDS (all genes) with caching -----------------------------
  dist_cache  <- file.path(cache_dir, "dist_all.RDS")
  mds_cache   <- file.path(cache_dir, "mds_all.RDS")
  message("Computing distances/MDS (all genes)…")
  if (!file.exists(dist_cache)) {
    dist_all <- parDist(t(qq_all), method = "euclidean", diag = FALSE, upper = FALSE, threads = threads)
    saveRDS(dist_all, dist_cache)
  } else dist_all <- readRDS(dist_cache)
  
  if (!file.exists(mds_cache)) {
    mds_all <- cmdscale(dist_all)
    colnames(mds_all) <- c("MDS1","MDS2")
    saveRDS(mds_all, mds_cache)
  } else mds_all <- readRDS(mds_cache)
  
  # ---- TPM (expressed-only) & QQ, Distances/MDS -----------------------------
  message("Filtering expressed genes & recomputing…")
  tpms_all2 <- tpms_all  # reuse
  is_expr <- rowSums(tpms_all2 >= tpm_threshold) >= (expr_prop * ncol(tpms_all2))
  qq_expr <- apply(tpms_all2[is_expr, , drop = FALSE], 2, function(x) qqnorm(x, plot.it = FALSE)$x)
  rownames(qq_expr) <- rownames(tpms_all2)[is_expr]
  assert_that(identical(colnames(qq_expr), meta_all$SAMPID))
  
  dist_cache_e  <- file.path(cache_dir, "dist_expressed.RDS")
  mds_cache_e   <- file.path(cache_dir, "mds_expressed.RDS")
  if (!file.exists(dist_cache_e)) {
    dist_expr <- parDist(t(qq_expr), method = "euclidean", diag = FALSE, upper = FALSE, threads = threads)
    saveRDS(dist_expr, dist_cache_e)
  } else dist_expr <- readRDS(dist_cache_e)
  
  if (!file.exists(mds_cache_e)) {
    mds_expr <- cmdscale(dist_expr)
    colnames(mds_expr) <- c("MDS1","MDS2")
    saveRDS(mds_expr, mds_cache_e)
  } else mds_expr <- readRDS(mds_cache_e)
  
  # ---- Plotting data frames --------------------------------------------------
  set.seed(seed)
  mds_df  <- merge(data.frame(SAMPID = rownames(mds_all), mds_all, row.names = NULL),
                   meta_all, by = "SAMPID")
  mds_e_df <- merge(data.frame(SAMPID = rownames(mds_expr), mds_expr, row.names = NULL),
                    meta_all, by = "SAMPID")
  
  # derive labels like your code
  mds_e_df$new_tissues2 <- ifelse(mds_df$SMTSD == "Brain - Frontal Cortex (BA9)", "GTEx PFC", "GTEx Brain")
  mds_e_df$new_tissues <- mds_df$SMTS
  mds_e_df$new_tissues <- ifelse(mds_e_df$new_tissues == "Brain", mds_e_df$new_tissues2, mds_e_df$new_tissues)
  mds_e_df$new_tissues <- gsub("PFC_msbb",   "MSBB PFC",   mds_e_df$new_tissues)
  mds_e_df$new_tissues <- gsub("PFC_lbp",    "LBP PFC",    mds_e_df$new_tissues)
  mds_e_df$new_tissues <- gsub("PFC_rosmap", "ROSMAP PFC", mds_e_df$new_tissues)
  
  mds_e_df$shape_group <- ifelse(mds_e_df$new_tissues %in% c("MSBB PFC","LBP PFC","ROSMAP PFC","GTEx PFC"),
                                 "PFC",
                          ifelse(mds_e_df$new_tissues == "GTEx Brain", "Non-PFC Brain", "Non-Brain"))
  shape_mapping <- c("PFC" = 17, "Non-PFC Brain" = 15, "Non-Brain" = 19)
  
  # shuffle colors but plot PFC on top (like your ordering trick)
  shuffled_levels <- sample(unique(mds_e_df$new_tissues))
  mds_e_df$new_tissues_shuffled <- factor(mds_e_df$new_tissues, levels = shuffled_levels)
  mds_e_df$shape_group <- factor(mds_e_df$shape_group, levels = c("Non-Brain","Non-PFC Brain","PFC"))
  # jittered plot (expressed genes)
  g4 <- ggplot(mds_e_df[order(mds_e_df$shape_group), ],
               aes(x = MDS1, y = MDS2, color = new_tissues_shuffled, shape = shape_group)) +
    geom_jitter(alpha = 0.5, width = 0.1, height = 0.1) +
    scale_color_discrete(name = "Tissue") +
    scale_shape_manual(name = "Tissue Type", values = shape_mapping) +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text  = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size = 2)),
           shape  = guide_legend(override.aes = list(size = 2)))
  
  # quick text and point plots (all genes), preserved from your script
  shuffled_idx <- sample(nrow(mds_df))
  g1 <- ggplot(mds_df[shuffled_idx, ], aes(x = MDS1, y = MDS2, color = SMTS, label = SMTS)) +
    geom_text(size = 1) + scale_color_discrete(name = "Tissue")
  g2 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = SMTSD)) +
    geom_point() + scale_color_discrete(name = "Tissue")
  g3 <- ggplot(subset(mds_df, SMTS %in% c("Brain","PFC_lbp","PFC_msbb","PFC_rosmap","Heart","Lung")),
               aes(x = MDS1, y = MDS2, color = SMTSD)) +
    geom_point() + scale_color_discrete(name = "Tissue")
  
  # ---- Save plots & table ----------------------------------------------------
  message("Saving plots and table…")
  
  # alt single (as in your script)
  pdf(pdf_out_alt, width = 10, height = 10)
  print(g4)
  dev.off()
  
  # export simplified table
  mds_export <- mds_e_df[, c("MDS1","MDS2","new_tissues_shuffled","shape_group")]
  colnames(mds_export) <- c("MDS1","MDS2","Tissue","Type")
  write_xlsx(mds_export, xlsx_out)
  
  # ---- Return a compact summary ---------------------------------------------
  invisible(list(
    n_samples = ncol(counts_all),
    n_genes_all = nrow(tpms_all),
    n_genes_expressed = sum(is_expr),
    pdf_alt  = pdf_out_alt,
    xlsx = xlsx_out,
    cache = list(dist_all = dist_cache, mds_all = mds_cache,
                 dist_expr = dist_cache_e, mds_expr = mds_cache_e)
  ))
}


res <- run_gtex_pfc_mds(
  # inputs
  rosmap_counts_path = "allcount_matrix_2023-04-11.txt",
  msbb_counts_path   = "allcount_matrix_2023-08-29.txt",
  lbp_counts_path    = "geneCounts_2022-07-26.txt",
  gtex_counts_gct_path       = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
  gtex_sample_attr_ds_path   = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
  gtex_subject_pheno_ds_path = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",

  # outputs
  output_dir = "/tmp/",
  cache_dir  = "/tmp",

  # params
  read_lengths = list(rosmap = 101, msbb = 100, lbp = 100, gtex = 76),
  tpm_threshold = 0.5,
  expr_prop = 0.10,
  threads = 30,
  seed = 42
)

