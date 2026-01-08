ml R/4.4.1 
R

####DE for bulkRNAseq and smri output = allSTRDEGresults_07.01.24.txt

##libraries 
rm(list=ls())
library(foreach)
library(doParallel)
registerDoParallel(cores=90)
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
wd = "/bulk_smri_noODC/"
expr = "geneCounts_2022-07-26.txt" #expression data on synapse
smri = "str_data_summarized_JAN172024.RDS" #see /data_processing/LBP/sMRI/smri_weights.R
cov = "cov_str_summarized_JAN172024.RDS" #meta info on synapse
meta = "info_all_withImaging_2022-07-26.txt" #meta info on synapse
phe = "lbpPheFile.txt" #see /misc/lbpPheFile.txt
results_file = "allSTRDEGresults_07.01.24.txt" #output
date="2024-06-23"
setwd(wd)

# ---------- Load RNA-seq counts and basic filtering ----------
genedata_org=fread(expr,data.table=FALSE)
genedata=as.data.frame(genedata_org[,2:ncol(genedata_org)])
rownames(genedata)=genedata_org[,1]
genedata=genedata[,colSums(genedata)>1.7e7] 

# ---------- Load sMRI-derived features and imaging covariates ----------
str_data_summarized=readRDS(smri)
str_data_summarized_weights=str_data_summarized$weights
str_data_summarized=str_data_summarized$E
cov_str_summarized=readRDS(cov)
cov_str_summarized=cov_str_summarized[match(colnames(str_data_summarized),cov_str_summarized$subject_id),]
assert_that(identical(cov_str_summarized$subject_id,colnames(str_data_summarized)))
cov_str_summarized$SurfaceHoles=str_data_summarized["SurfaceHoles",]
cov_str_summarized$rhSurfaceHoles=str_data_summarized["rhSurfaceHoles",]
cov_str_summarized$lhSurfaceHoles=str_data_summarized["lhSurfaceHoles",]

# ---------- Load sample metadata and phenotype, align to RNA-seq ----------
info_all=fread(meta,sep="\t",header=TRUE,data.table=FALSE)
pheTable=fread(phe,data.table=FALSE)
info_all=merge(info_all,pheTable,by.x="IID_ISMMS",by.y="iid",suffixes=c("","All"))
assert_that(identical(info_all$SAMPLE_ISMMS[match(colnames(genedata),info_all$SAMPLE_ISMMS)],colnames(genedata)))
info_all2=info_all[match(colnames(genedata),info_all$SAMPLE_ISMMS),]
info_all2=as.data.frame(info_all2)

# ---------- Sanity checks for imaging↔RNA subject overlap ----------
setdiff(colnames(str_data_summarized),info_all2$IID_ISMMS)
setdiff(info_all2$IID_ISMMS,colnames(str_data_summarized))

# ---------- Align everything to a common set and bind imaging features ----------
IID_withRNAandImaging=intersect(colnames(str_data_summarized),info_all2$IID_ISMMS)
fullInfo=merge(cov_str_summarized,info_all2,by.x="subject_id",by.y="IID_ISMMS",all.y=TRUE,suffixes=c(".imaging",".RNA"))
assert_that(identical(fullInfo$SAMPLE_ISMMS[match(colnames(genedata),fullInfo$SAMPLE_ISMMS)],colnames(genedata)))
fullInfo2=fullInfo[match(colnames(genedata),fullInfo$SAMPLE_ISMMS),]
fullInfo2=as.data.frame(fullInfo2)
fullInfo2=fullInfo2[fullInfo2$subject_id %in% colnames(str_data_summarized),]
genedata=genedata[,fullInfo2$SAMPLE_ISMMS]
assert_that(identical(colnames(genedata),fullInfo2$SAMPLE_ISMMS ))
str_data_summarized=str_data_summarized[,fullInfo2$subject_id]
assert_that(identical(colnames(str_data_summarized),fullInfo2$subject_id))
fullInfo2=cbind(fullInfo2,t(str_data_summarized))
assert_that(identical(fullInfo2$SAMPLE_ISMMS,colnames(genedata)))

# ---------- Filter lowly-expressed genes (CPM ≥1 in ≥10% of samples) ----------
isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
table(isexpr)

# ---------- Build DGEList and normalize ----------
genesAll <- DGEList(counts=genedata[isexpr,])#, group=info$disease) 
genesAll <- calcNormFactors(genesAll) 

dim(genesAll)
vobj=voom(genesAll,plot=F)

# ---------- Use UTY/XIST to infer expressed sex ----------
UTY="ENSG00000183878"
Xist="ENSG00000229807"

UTY_index=grep(UTY,rownames(vobj))
Xist_index=grep(Xist,rownames(vobj))

fullInfo2$expressedSex="male"
fullInfo2$expressedSex[vobj$E[UTY_index,]<vobj$E[Xist_index,]]="female"
fullInfo2$expressedSex[is.na(fullInfo2$SAMPLE_ISMMS)]=NA
rownames(fullInfo2)=fullInfo2$SAMPLE_ISMMS

# ---------- Scale some covariates/features for stability ----------
fullInfo2$EstimatedTotalIntraCranialVol_scaled=scale(fullInfo2$EstimatedTotalIntraCranialVol)
fullInfo2$SurfaceHoles_scaled=scale(fullInfo2$SurfaceHoles)
fullInfo2$lh_rostralmiddlefrontal_thickness_scaled=scale(fullInfo2$lh_rostralmiddlefrontal_thickness)
fullInfo2$rh_caudalanteriorcingulate_thickness_scaled=scale(fullInfo2$rh_caudalanteriorcingulate_thickness)

# ---------- Parallelize ----------
ncpus <- 36
nprocs <- 36                                          
threads_per_proc = floor(ncpus / nprocs)
Sys.setenv(OMP_NUM_THREADS=as.character(threads_per_proc))

# ---------- voomWithDreamWeights design  ----------
formVoom = formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + (1|expressedSex) + scale(RNASeqMetrics_PCT_INTRONIC_BASES) + (1|mymet_tissue) + SurfaceHoles"))

# Construct a file-safe name from the formula for voomWithDreamWeights result
forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
if(length(which(paste0(forFileNameVoom)==list.files("expression/")))==0){
  system(paste0("mkdir expression/",forFileNameVoom))
  vobjDream = voomWithDreamWeights( genesAll, formVoom, fullInfo2, save.plot=TRUE, BPPARAM = MulticoreParam(nprocs))
  saveRDS(vobjDream,file=paste0("expression/",forFileNameVoom,"/voomWithDreamWeights.RDS"))
}else{
  vobjDream=readRDS(paste0("expression/",forFileNameVoom,"/voomWithDreamWeights.RDS"))
}
assert_that(identical(fullInfo2$SAMPLE_ISMMS,colnames(vobjDream)))

# ---------- Per-feature DREAM fits using imaging features as covariates ----------
# Loop over all imaging features except those containing "SurfaceHoles" (already in model)
res_final=foreach(feature=grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE))%do%{
	cat(feature,"\n")
  if(feature=="EstimatedTotalIntraCranialVol"){
    feature="EstimatedTotalIntraCranialVol_scaled"
  }
  # if(feature=="lh_rostralmiddlefrontal_thickness" | feature=="rh_caudalanteriorcingulate_thickness" | feature=="rh_transversetemporal_thickness" | feature=="Left_Putamen" | feature=="Left_Pallidum" | feature=="Right_Pallidum"){
  #   feature=paste0(feature,"_scaled")
  # }
  cat(feature," renamed\n")
	form=formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + (1|expressedSex) + RNASeqMetrics_PCT_INTRONIC_BASES + (1|mymet_tissue) + SurfaceHoles + ",feature))
  if(feature=="lh_rostralmiddlefrontal_thickness_scaled" | feature=="Left_vessel"){
    form=formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + (1|expressedSex) + RNASeqMetrics_PCT_INTRONIC_BASES + (1|mymet_tissue) + SurfaceHoles_scaled + ",feature))
  }
	design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(form)[2],fixed=TRUE))),fullInfo2)
	vobjDreamSubset=vobjDream[,match(rownames(design),colnames(vobjDream))]
	fullInfo2Subset=fullInfo2[match(rownames(design),fullInfo2$SAMPLE_ISMMS),]
	assert_that(identical(fullInfo2Subset$SAMPLE_ISMMS,colnames(vobjDreamSubset)))

	forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(form)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

	if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
		fitDream=tryCatch({dream(vobjDreamSubset, form, fullInfo2Subset, BPPARAM = MyBatchtoolsLSFParam(conf.file = "/.batchtools.conf.R", workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = 15000)))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror while running dream on ",feature,": "),e); NULL}) #, ddf="Kenward-Roger")
		if(is.null(fitDream)==FALSE){
		              system(paste0("mkdir DE/",forFileName))
		              saveRDS(fitDream,file=paste0("DE/",forFileName,"/fitDream.RDS"))
		}
	}
}

# ---------- Read back per-feature dream fits and collect topTable results ----------
res_final=foreach(feature=grep("SurfaceHoles",rownames(str_data_summarized),value=TRUE,invert=TRUE))%do%{
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
  form=formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + (1|expressedSex) + RNASeqMetrics_PCT_INTRONIC_BASES + (1|mymet_tissue) + SurfaceHoles + ",feature))
  if(feature=="lh_rostralmiddlefrontal_thickness_scaled" | feature=="Left_vessel"){
    form=formula(paste0("~ (1|subject_id) + (1|mymet_depletionbatch) + mymet_rin + (1|mymet_bank) + RNASeqMetrics_PCT_INTRONIC_BASES + (1|mymet_tissue) + SurfaceHoles_scaled + ",feature))
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




