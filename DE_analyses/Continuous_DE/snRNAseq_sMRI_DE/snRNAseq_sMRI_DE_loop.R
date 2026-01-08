ml R/4.4.1
R

##code for processing and running snRNAseq sMRI DE - output top_table_results_sc_DE_results_1.15.25.RDS and pi_sc_DE_results_1.15.25.RDS

# ---------- packages ----------
library(dreamlet)
library(assertthat)
library(foreach)
library(doParallel)
library(Seurat)
library(BiocParallel)
register(MulticoreParam(workers = 30))
library(limma)
library(data.table)

#  [1] Seurat_5.1.0                SeuratObject_5.0.2
#  [3] sp_2.1-4                    doParallel_1.0.17
#  [5] iterators_1.0.14            foreach_1.5.2
#  [7] assertthat_0.2.1            dreamlet_1.4.1
#  [9] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [11] Biobase_2.66.0              GenomicRanges_1.56.0
# [13] GenomeInfoDb_1.42.1         IRanges_2.40.1
# [15] S4Vectors_0.44.0            BiocGenerics_0.52.0
# [17] MatrixGenerics_1.18.0       matrixStats_1.4.1
# [19] variancePartition_1.36.3    BiocParallel_1.40.0
# [21] limma_3.62.1                ggplot2_3.5.1


# ---------- Project paths and inputs ----------
setwd("/add_sex/")

str_data_summarized = "str_data_summarized_JAN172024.RDS" #see /data_processing/LBP/sMRI/smri_weights.R
lbp = "cpm_log2_seurat_scaled-nCount-nFeature-2000g-pca15-harmony-umap-clust0p6-anno.RDS" #see syanpse for expression 
cov = "Liv_samples_metadata_for_Anina.RDS" # see synapse for meta
folder = "/add_sex/"  #folder
cov_str_summarized = "cov_str_summarized_JAN172024.RDS" #see synapse for meta
results = "/add_sex/" #output for results
toptable_results = "top_table_results_sc_DE_results_1.15.25.RDS" #output
pi1_results = "pi_sc_DE_results_1.15.25.RDS" #output
DE_output = "/DE/*RDS" #output

# ---------- Load objects ----------
lbp = readRDS(lbp) 
cov = readRDS(cov) 
living_ids = cov$ids

# ---------- Convert to SingleCellExperiment and set assay ----------
pb <- list()
sce <- as.SingleCellExperiment(lbp)
assayNames(sce)[1] = "counts"

# ---------- Pseudobulk aggregation ----------
file = paste0(folder, "expression/lbp_01.14.25_pseudobulk.RDS")
if (!file.exists(file)){
    pb = aggregateToPseudoBulk(sce,
        assay = "counts", 
        cluster_id = "CellTypeAWC",
        sample_id = "Sample_ID",
        BPPARAM = MulticoreParam(30))
    saveRDS(pb, file)
}else{
    pb=readRDS(file)
}

#check
assayNames(pb)

# ---------- Restrict to living cohort ----------
subset_liv <- pb[, c("PT-0236R", "PT-0237R", "PT-0256L", "PT-0257R", "PT-0261L", 
"PT-0262L", "PT-0264L", "PT-0264R", "PT-0267L", "PT-0267R", "PT-0280L", 
"PT-0280R", "PT-0282L", "PT-0282R", "PT-0283L", "PT-0284R", "PT-0285L", 
"PT-0285R", "PT-0286L", "PT-0286R", "PT-0287L", "PT-0288L", "PT-0289L", 
"PT-0289R", "PT-0290L", "PT-0290R", "PT-0292L", "PT-0292R", "PT-0295L", 
"PT-0296L", "PT-0297L")]

# ---------- Derive total nCells per sample and attach to colData ----------
ncells_variable=as.data.frame((int_colData(subset_liv)$n_cells))
ncells_variable=t(ncells_variable)

# Calculate row-wise sums and create a new row
sums_row <- rowSums(ncells_variable)
sums_row = as.data.frame(sums_row)
rownames(sums_row) <- gsub("\\.", "-", rownames(sums_row))
colnames(sums_row)="ncells"

## Merge total nCells into colData for modeling as a covariate
colData(subset_liv)$ncells = sums_row

# ---------- Bring in imaging-derived features and align with samples ----------
str_data_summarized=readRDS(str_data_summarized)
str_data_summarized_weights=str_data_summarized$weights
str_data_summarized=str_data_summarized$E
cov_str_summarized=readRDS(cov_str_summarized)
cov_str_summarized=cov_str_summarized[match(colnames(str_data_summarized),cov_str_summarized$subject_id),]
assert_that(identical(cov_str_summarized$subject_id,colnames(str_data_summarized)))
cov_str_summarized$SurfaceHoles=str_data_summarized["SurfaceHoles",]
cov_str_summarized$rhSurfaceHoles=str_data_summarized["rhSurfaceHoles",]
cov_str_summarized$lhSurfaceHoles=str_data_summarized["lhSurfaceHoles",]

# ---------- Build unified metadata (single-cell + imaging) ----------
info_all2=DataFrame(subset_liv@colData)
info_all2$Sample_ID=rownames(info_all2)

subset_liv_order=colnames(subset_liv@assays@data$Exc2)

# Merge imaging covariates (keyed by subject_id) with single-cell pseudobulk metadata (keyed by Indv_ID)
fullInfo=merge(cov_str_summarized,info_all2,by.x="subject_id",by.y="Indv_ID",all.y=TRUE,suffixes=c(".imaging",".single_cell"))

## Sort rows to match order of pseudobulk assay columns
fullInfo2=fullInfo[match(colnames(subset_liv@assays@data$Exc),fullInfo$Sample_ID),]

# Reduce imaging matrix to the aligned subjects and double-check ordering
str_data_summarized=str_data_summarized[,fullInfo2$subject_id]
assert_that(identical(colnames(str_data_summarized),fullInfo2$subject_id))
#TRUE

fullInfo2=cbind(fullInfo2,t(str_data_summarized))
assert_that(identical(fullInfo2$Sample_ID,colnames(subset_liv@assays@data$Exc)))
#[1] TRUE

# ---------- Scale some covariates/features for stability ----------
rownames(fullInfo2)=fullInfo2$Sample_ID

fullInfo2$EstimatedTotalIntraCranialVol_scaled=scale(fullInfo2$EstimatedTotalIntraCranialVol)
fullInfo2$SurfaceHoles_scaled=scale(fullInfo2$SurfaceHoles)
fullInfo2$lh_rostralmiddlefrontal_thickness_scaled=scale(fullInfo2$lh_rostralmiddlefrontal_thickness)
fullInfo2$rh_caudalanteriorcingulate_thickness_scaled=scale(fullInfo2$rh_caudalanteriorcingulate_thickness)

#### ---------- Add enriched metadata back into subset_liv colData ----------
mergedMeta=dplyr::inner_join(tibble::rownames_to_column(data.frame(subset_liv@colData)),tibble::rownames_to_column(data.frame(fullInfo2)))
assert_that(identical(mergedMeta$rowname,rownames(data.frame(subset_liv@colData))))
#true

rownames(mergedMeta)=mergedMeta$rowname
mergedMeta$rowname=NULL
subset_liv@colData=DataFrame(mergedMeta)

# saveRDS(subset_liv, "pre_normalize_subset_liv_1.24.24.RDS")

# ---------- Modeling formula for normalization (dreamlet::processAssays) ----------
form <-  formula(paste0(" ~ (1|Batch) + (1|mymet_PD) + (1|Indv_ID) + scale(Mean_Reads_per_Cell) + scale(Fraction_Reads_in_Cells) + scale(ncells) + scale(SurfaceHoles)"))

# ---------- Normalize per assay ----------
file = paste0(folder, "expression/lbp_01.14.25_pseudobulk_normalized.RDS")
if (!file.exists(file)){
    res.proc = processAssays(subset_liv, form, BPPARAM = MulticoreParam(20))
    saveRDS(res.proc, file)
}else{
    res.proc=readRDS(file)
}

####add fullInfo2 back into subset_liv
mergedMeta=dplyr::inner_join(tibble::rownames_to_column(data.frame(subset_liv@colData)),tibble::rownames_to_column(data.frame(fullInfo2)))
assert_that(identical(mergedMeta$rowname,rownames(data.frame(subset_liv@colData))))

#true
rownames(mergedMeta)=mergedMeta$rowname
mergedMeta$rowname=NULL
subset_liv@colData=DataFrame(mergedMeta)


# ---------- Per-feature DREAM fits using imaging features as covariates ----------
res_final=foreach(feature=grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE))%do%{
   cat(feature,"\n")
  if(feature=="EstimatedTotalIntraCranialVol"){
    feature="EstimatedTotalIntraCranialVol_scaled"
    cat(feature," renamed\n")
  }
  
  form=formula(paste0("~(1|Batch) + (1|mymet_PD) + (1|Indv_ID) + scale(Mean_Reads_per_Cell) + scale(Fraction_Reads_in_Cells) + scale(ncells) + scale(SurfaceHoles) + ", feature))

  if(feature=="lh_rostralmiddlefrontal_thickness" | feature=="Left_vessel"){
    form=formula(paste0("~ (1|Batch) + (1|mymet_PD) + (1|Indv_ID) + scale(Mean_Reads_per_Cell) + scale(Fraction_Reads_in_Cells) + scale(ncells) + scale(SurfaceHoles) + scale(",feature, ")"))
  }

  forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(form)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
    # file = paste0(folder, "DE/", forFileName)
    # if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
    folder = results
    file = paste0(folder, "DE/", forFileName, ".RDS")
    if (!file.exists(file)){
    fit = dreamlet(res.proc, form, BPPARAM = MulticoreParam(30))
    saveRDS(fit, file)
    }else{
    fit=readRDS(file)
    }
}


# ---------- Collect all per-feature fit files ----------
flist <- Sys.glob(DE_output)

# ---------- Summarize DE: pi1 and # of FDR-significant genes per assay ----------
myres <- c()
# myres_for_prot <- c()
failures <- c()
for (i in flist){
  name <- gsub(".RDS", "", unlist(strsplit(basename(i), split="SurfaceHoles_")))[2]
  data <- readRDS(i)
  if (!is.null(coefNames(data))){
    data <- topTable(data, coef=as.character(coefNames(data)[6]), num=Inf)
    data <- data.table(feature=name, as.data.table(data))
    for (j in unique(data$assay)){
      cur <- data[assay==j]
      p1 <- 1 - propTrueNull(cur$P.Value)
      nde <- nrow(cur[adj.P.Val<=0.05])
      add <- data.table( feature=name, assay=j, pi1=p1, ndeg=nde)
      myres <- rbind(myres, add)
      print(add, row.names=FALSE, col.names="none")
      print("", row.names=FALSE)
    }
  } else {
    failures <- c(failures, i)
  }
}

saveRDS(myres, pi1_results)


# ---------- Save concatenated topTables across features/assays ----------
myres2 <- c()
# myres2_for_prot <- c()
failures <- c()
for (i in flist){
  name <- gsub(".RDS", "", unlist(strsplit(basename(i), split="SurfaceHoles_")))[2]
  data <- readRDS(i)
  if (!is.null(coefNames(data))){
    print(name)
    data <- topTable(data, coef=as.character(coefNames(data)[6]), num=Inf)
    data <- data.table(feature=name, as.data.table(data))
    myres2 = rbind(data, myres2)

    }}

saveRDS(myres2, toptable_results)
