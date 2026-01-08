rm(list=ls())
##################################################
#load libraries
##################################################


library(geneLenDataBase)
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(gridExtra)
library(PRROC)
library(BiocParallel)
library(Matrix)
library(lme4)
register(SnowParam(45))
library(assertthat)
library(matrixStats)

#WEIGHTS
MPRAGE_cov=cov_str[cov_str$SeriesDescription2=="MPRAGE",]
MPRAGE_contrast_cov=cov_str[cov_str$SeriesDescription2=="MPRAGE_contrast",]
FSPGR_cov=cov_str[cov_str$SeriesDescription2=="FSPGR",]
FSPGR_contrast_cov=cov_str[cov_str$SeriesDescription2=="FSPGR_contrast",]

MPRAGE=str_data_norm[,match(MPRAGE_cov$scan_id,colnames(str_data_norm))]
assert_that(identical(as.character(MPRAGE_cov$scan_id),colnames(MPRAGE)))
colnames(MPRAGE)=MPRAGE_cov$subject_id


MPRAGE_contrast=str_data_norm[,match(MPRAGE_contrast_cov$scan_id,colnames(str_data_norm))]
assert_that(identical(as.character(MPRAGE_contrast_cov$scan_id),colnames(MPRAGE_contrast)))
colnames(MPRAGE_contrast)=MPRAGE_contrast_cov$subject_id

FSPGR=str_data_norm[,match(FSPGR_cov$scan_id,colnames(str_data_norm))]
assert_that(identical(as.character(FSPGR_cov$scan_id),colnames(FSPGR)))
colnames(FSPGR)=FSPGR_cov$subject_id

FSPGR_contrast=str_data_norm[,match(FSPGR_contrast_cov$scan_id,colnames(str_data_norm))]
assert_that(identical(as.character(FSPGR_contrast_cov$scan_id),colnames(FSPGR_contrast)))
colnames(FSPGR_contrast)=FSPGR_contrast_cov$subject_id

MPRAGE_arrayWeight=arrayWeights(object=MPRAGE)

#remove no variance
# Assuming str_data_norm is your data frame
rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities", "Left_non_WM_hypointensities", "Right_non_WM_hypointensities")
rows_indices_to_remove <- which(rownames(str_data_norm) %in% rows_to_remove)
rows_to_remove <- c(35, 36, 38, 39)
str_data_norm2 <- str_data_norm[-rows_to_remove, ]
###

varsMPRAGE=rowVars(as.matrix(MPRAGE))
varsMPRAGE_contrast=rowVars(as.matrix(MPRAGE_contrast))
varsFSPGR=rowVars(as.matrix(FSPGR))
varsFSPGR_contrast=rowVars(as.matrix(FSPGR_contrast))

# do the weights after adjusting for batch effects using comBat
imageTypeWeights=c(MPRAGE=1/median(varsMPRAGE),MPRAGE_contrast=1/median(varsMPRAGE_contrast),FSPGR=1/median(varsFSPGR),FSPGR_contrast=1/median(varsFSPGR_contrast))
imageTypeWeights=exp(scale(log(imageTypeWeights),center=TRUE,scale=FALSE))[,1]

#make sure scans are on the same day for the ones we do weighted mean for
booleanContribution=matrix(data=FALSE,ncol=4,nrow=length(unique(cov_str$subject_id)),dimnames=list(unique(cov_str$subject_id),names(imageTypeWeights)))
tmp=c()
for(subject in unique(cov_str$subject_id)){
	cols=cov_str$subject_id==subject & cov_str$SeriesDescription2 %in% names(imageTypeWeights)
	tmp=c(tmp,paste0(names(imageTypeWeights[cov_str$SeriesDescription2[cols]]),collapse="_"))
	booleanContribution[subject,names(imageTypeWeights[cov_str$SeriesDescription2[cols]])]=TRUE
}

fractionContribution=scale(booleanContribution,scale=1/imageTypeWeights,center=FALSE)
fractionContribution=fractionContribution/rowSums(fractionContribution)


selected_samples=cov_str$scan_id[cov_str$SeriesDescription2 %in% names(imageTypeWeights)]
selected_samples=as.character(selected_samples)
str_data_norm2=str_data_norm2[,selected_samples]
cov_str=cov_str[match(selected_samples,cov_str$scan_id),]
assert_that(identical(as.character(cov_str$scan_id),colnames(str_data_norm2)))
#[1] TRUE

imageTypeWeightsMatrix=expandAsMatrix(imageTypeWeights[cov_str$SeriesDescription2],dim=dim(str_data_norm2))
to_remove <- c(as.numeric(names(which(table(which(is.na(cov_str),arr.ind=TRUE)[,2])==nrow(cov_str)))),which(apply(cov_str, 2,function(x){length(unique(x[is.na(x)==FALSE]))==1})))
if(length(to_remove)>0){
cov_str2<- cov_str[,-to_remove]}
rownames(cov_str)=cov_str$scan_id

###IMPORTANT - remove effect of seriesdescription while keep subject effect and use weight o frach technology 
fitList = fitVarPartModel(str_data_norm2, ~ (1|subject_id) + (1|SeriesDescription2), cov_str, weightsMatrix=imageTypeWeightsMatrix, fxn = function (fit) {
            ## Subtract batch effect, keep sample effect & residuals
            y_batch_corrected <- variancePartition::get_prediction( fit, ~ 1 + (1|subject_id)) + resid(fit)
            subject_id <- model.frame(fit)$subject_id
            ## Get the weights used in the fit
            w <- model.frame(fit)[["(weights)"]]
            # use lm to get the mean and standard error for each sample
            fit2 <- lm(
                y_batch_corrected ~ 0 + subject_id,
                weights = w
            )
            predictFit2 <- predict.lm(fit2, se.fit=TRUE)
            ## Coefs are summarized expression values
            fit2_coef <- coef(fit2)
            ## Coef inverse variances are summarised weights
            fit2_coef_var <- exp(scale(log(diag(vcov(fit2))),center=TRUE,scale=FALSE))
            ## Strictly speaking this should be a data frame, but we
            ## use a matrix because it's more memory-efficient
            res <- cbind(
                expr = fit2_coef,
                weight = 1/fit2_coef_var
            )
            rownames(res) <- fit2$xlevels$subject_id
            res
        }, 
    showWarnings=FALSE, quiet=TRUE, BPPARAM=MulticoreParam())

values=t(sapply(fitList,function(x){x[,1]}))
weights=t(sapply(fitList,function(x){x[,2]}))

str_data_summarized <- new("EList")
str_data_summarized$E=values
str_data_summarized$weights=weights

saveRDS(str_data_summarized,"str_data_summarized_JAN172024.RDS")
