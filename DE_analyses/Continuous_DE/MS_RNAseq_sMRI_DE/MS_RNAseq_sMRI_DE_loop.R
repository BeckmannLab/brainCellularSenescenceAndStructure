ml R/4.4.1 
R

####DE for bulkRNAseq and smri output = allSTRDEGresults_01.18.24.txt

##libraries 
rm(list=ls())
library(foreach)
library(doParallel)
registerDoParallel(cores=30)
library(variancePartition)
library(memoise)
library(dplyr)
library(data.table)
library(limma)
library(edgeR)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(fs)

cachedir <- "/directory" # Change this to whatever you're using
dir_create(cachedir)
file_chmod(cachedir, "go-rwx")
fc <- cache_filesystem(cachedir)

## Memoise some functions using the above cache
voom <- memoise(limma::voom, cache = fc)
voomWithDreamWeights <- memoise(
    variancePartition::voomWithDreamWeights,
    cache = fc,
    omit_args = "BPPARAM")
dream <- memoise(
    variancePartition::dream,
    cache = fc,
    omit_args = "BPPARAM")

MyBatchtoolsLSFParam <- function(..., registryargs = batchtoolsRegistryargs(),
                                  conf.file = "/.batchtools.conf.R") {
    if (missing(conf.file)) {
        return(BatchtoolsParam(..., registryargs = registryargs))
    }
    registryargs$conf.file <- conf.file
    bt_conf <- within(list(), {
        source(conf.file, local = TRUE)
    })
    myargs <- list(
        cluster = "lsf",
        registryargs = registryargs,
        ...
    )
    param <- BatchtoolsParam(
        cluster = "lsf",
        registryargs = registryargs,
        ...
    )
    if (!is.null(bt_conf$default.resources)) {
        merged_resources <- bt_conf$default.resources
        merged_resources[names(param$resources)] <- param$resources
        param$resources <- merged_resources
    }
    bpstart(param)
    if (!is.null(bt_conf$cluster.functions)) {
        param$registry$cluster.functions <- bt_conf$cluster.functions
    }
    param
}

#paths EDIT
date="2024-01-17"
setwd("/proteomics/")
str_data_summarized = "str_data_summarized_JAN172024.RDS" #see /data_processing/LBP/sMRI/smri_weights.R
protdata_org = "expression_liv_01.16.24.RDS" #see synapse protein 
cov_str_summarized = "cov_str_summarized_JAN172024.RDS" #meta info on synapse
meta = "covariates_liv_01.16.24.RDS"  #meta info on synapse
phe = "lbpPheFile.txt" #see /misc/lbpPheFile.txt
results_file = "allSTRDEGresults_01.18.24.txt"
# ---- Load protein data and drop outlier sample ----
protdata_org=readRDS(protdata_org) 
protdata=protdata_org[,-which(colnames(protdata_org)=="LBPSEMA4BRAIN618")] ##remove outlier 

## ---- Imaging set up and formatting ----
str_data_summarized=readRDS(str_data_summarized)
str_data_summarized_weights=str_data_summarized$weights
str_data_summarized=str_data_summarized$E
cov_str_summarized=readRDS(cov_str_summarized)
cov_str_summarized=cov_str_summarized[match(colnames(str_data_summarized),cov_str_summarized$subject_id),]
assert_that(identical(cov_str_summarized$subject_id,colnames(str_data_summarized)))
#TRUE
cov_str_summarized$SurfaceHoles=str_data_summarized["SurfaceHoles",]
cov_str_summarized$rhSurfaceHoles=str_data_summarized["rhSurfaceHoles",]
cov_str_summarized$lhSurfaceHoles=str_data_summarized["lhSurfaceHoles",]

# ---- Merge imaging covariates with phenotype/covariates, align with proteomics ----
info_all=readRDS(meta)
pheTable=fread(phe,data.table=FALSE)
info_all=merge(info_all,pheTable,by.x="IID_ISMMS",by.y="iid",suffixes=c("","All"))
assert_that(identical(info_all$SAMPLE_ISMMS[match(colnames(protdata),info_all$SAMPLE_ISMMS)],colnames(protdata)))
#TRUE

info_all2=info_all[match(colnames(protdata),info_all$SAMPLE_ISMMS),]
info_all2=as.data.frame(info_all2)

dput(setdiff(colnames(str_data_summarized),info_all2$IID_ISMMS))

dput(setdiff(info_all2$IID_ISMMS,colnames(str_data_summarized)))

IID_withRNAandImaging=intersect(colnames(str_data_summarized),info_all2$IID_ISMMS)
##155 IDs in common
fullInfo=merge(cov_str_summarized,info_all2,by.x="subject_id",by.y="IID_ISMMS",all.y=TRUE,suffixes=c(".imaging",".proteomics"))
assert_that(identical(fullInfo$SAMPLE_ISMMS[match(colnames(protdata),fullInfo$SAMPLE_ISMMS)],colnames(protdata)))
#[1] TRUE
fullInfo2=fullInfo[match(colnames(protdata),fullInfo$SAMPLE_ISMMS),]
fullInfo2=as.data.frame(fullInfo2)

## Align to subjects present in imaging (restrict to intersection)
fullInfo2=fullInfo2[fullInfo2$subject_id %in% colnames(str_data_summarized),]
protdata=protdata[,fullInfo2$SAMPLE_ISMMS]
assert_that(identical(colnames(protdata),fullInfo2$SAMPLE_ISMMS ))

str_data_summarized=str_data_summarized[,fullInfo2$subject_id]
assert_that(identical(colnames(str_data_summarized),fullInfo2$subject_id))
#[1] TRUE

# Bind imaging metrics (as columns) to covariate frame; re-align protein matrix
fullInfo2=cbind(fullInfo2,t(str_data_summarized))
protdata2=protdata[,match(fullInfo2$SAMPLE_ISMMS,colnames(protdata))]

assert_that(identical(fullInfo2$SAMPLE_ISMMS,colnames(protdata2)))
#[1] TRUE
ProtAll =protdata2

# ---------- Scale some covariates/features for stability ----------
rownames(fullInfo2)=fullInfo2$SAMPLE_ISMMS
fullInfo2$EstimatedTotalIntraCranialVol_scaled=scale(fullInfo2$EstimatedTotalIntraCranialVol)
fullInfo2$SurfaceHoles_scaled=scale(fullInfo2$SurfaceHoles)
fullInfo2$lh_rostralmiddlefrontal_thickness_scaled=scale(fullInfo2$lh_rostralmiddlefrontal_thickness)
fullInfo2$rh_caudalanteriorcingulate_thickness_scaled=scale(fullInfo2$rh_caudalanteriorcingulate_thickness)

# ---------- Parallelize ----------
ncpus <- 36
nprocs <- 36                                          
threads_per_proc = floor(ncpus / nprocs)
Sys.setenv(OMP_NUM_THREADS=as.character(threads_per_proc))
register(SerialParam())

# ---------- voomWithDreamWeights design  ----------
# formVoom = formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + (1|expressedSex) + RNASeqMetrics_PCT_INTRONIC_BASES + ODC + (1|mymet_tissue) + SurfaceHoles"))
formVoom = formula(paste0("~  (1|Bank) + (1|SINAI_SizeGuesstimate) + (1|subject_id) + SurfaceHoles"))

# Construct a file-safe name from the formula for voomWithDreamWeights result
forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

# ---------- Per-feature DREAM fits using imaging features as covariates ----------
res_final=foreach(feature=grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE))%dopar%{
   cat(feature,"\n")
  if(feature=="EstimatedTotalIntraCranialVol"){
    feature="EstimatedTotalIntraCranialVol_scaled"
  }
  cat(feature," renamed\n")
  form=formula(paste0("~(1|Bank) + (1|SINAI_SizeGuesstimate) + (1|subject_id) + SurfaceHoles + ",feature))

  if(feature=="lh_rostralmiddlefrontal_thickness_scaled" | feature=="Left_vessel"){
    form=formula(paste0("~  (1|Bank) + (1|SINAI_SizeGuesstimate) + (1|subject_id) + SurfaceHoles_scaled + ",feature))
  }

  design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(form)[2],fixed=TRUE))),fullInfo2)

  ProtAllSubset=ProtAll[,match(rownames(design),colnames(ProtAll))]
  fullInfo2Subset=fullInfo2[match(rownames(design),fullInfo2$SAMPLE_ISMMS),]
  assert_that(identical(fullInfo2Subset$SAMPLE_ISMMS,colnames(ProtAllSubset)))

  forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(form)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

  if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
    fitDream=tryCatch({dream(ProtAllSubset, form, fullInfo2Subset, BPPARAM = MyBatchtoolsLSFParam(conf.file = "/.batchtools.conf.R", workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = "15GB")))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror while running dream on ",feature,": "),e); NULL}) #, ddf="Kenward-Roger")
    if(is.null(fitDream)==FALSE){
                  system(paste0("mkdir DE/",forFileName))
                  saveRDS(fitDream,file=paste0("DE/",forFileName,"/fitDream.RDS"))
    }
  }
}

# ---------- Read back per-feature dream fits and collect topTable results ----------
res_final=foreach(feature=grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE))%dopar%{
# for(feature in grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE)){
  cat(feature,"\n")
  if(feature=="EstimatedTotalIntraCranialVol"){
    feature="EstimatedTotalIntraCranialVol_scaled"
    cat(feature," renamed\n")
  }
  # if(feature=="lh_rostralmiddlefrontal_thickness"){
  #   feature="lh_rostralmiddlefrontal_thickness_scaled"
  #   cat(feature," renamed\n")
  # }
  form=formula(paste0("~ (1|Bank) + (1|SINAI_SizeGuesstimate) + (1|subject_id) + SurfaceHoles + ",feature))
  if(feature=="lh_rostralmiddlefrontal_thickness_scaled" | feature=="Left_vessel"){
    form=formula(paste0("~ (1|Bank) + (1|SINAI_SizeGuesstimate) + (1|subject_id) + SurfaceHoles_scaled + ",feature))
  }
  forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(form)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
  fitDream=readRDS(paste0("DE/",forFileName,"/fitDream.RDS"))
    res1 = topTable(fitDream, coef=feature, number=Inf, sort.by="P")
    res1$feature=feature
    res1$geneID=rownames(res1)
    res1
}


# Bind all topTables together
resFinal=do.call(rbind,res_final)
fwrite(resFinal,file=results_file,col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

