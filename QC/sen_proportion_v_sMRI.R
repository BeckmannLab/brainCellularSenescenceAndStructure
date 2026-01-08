#setup 
rm(list=ls())
library(data.table)
library(ggplot2)
map <- fread("~/gene_ids_ensembl2symbol_fromHUGO_10JUN2020.tsv", na=c("","NA"), col.names=c("symbol", "gene"))
map <- map[!is.na(gene)]

##data
scores <- list(
ast = #aucell scores
exc = #aucell scores,
int = #aucell scores,
mg = #aucell scores,
nonneu = #aucell scores,
oli = #aucell scores,
opc = #aucell scores,
all = readRDS("celltype_146_senescence_status.RDS"))
str_data_summarized=readRDS("/str_data_summarized_JAN172024.RDS")
str_data_summarized_weights=str_data_summarized$weights
str_data_summarized=str_data_summarized$E
cov_str_summarized=readRDS("cov_str_summarized_JAN172024.RDS")
cov_str_summarized=cov_str_summarized[match(colnames(str_data_summarized),cov_str_summarized$subject_id),]
cov_str_summarized$SurfaceHoles=str_data_summarized["SurfaceHoles",]
cov_str_summarized$rhSurfaceHoles=str_data_summarized["rhSurfaceHoles",]
cov_str_summarized$lhSurfaceHoles=str_data_summarized["lhSurfaceHoles",]

##collapse scores
myoutput <- c()
for (i in names(scores)){
  cur <- as.data.table(scores[[i]], keep.rownames="cell")
  cur[,sample:=tstrsplit(cell, split="_", keep=1L)]
  ncells <- cur[,list(ncell=.N),sample] 
  ntruth <- cur[sen==TRUE,list(ntrue=.N),sample] 
  mer <- merge(ncells, ntruth, all=TRUE)
  mer[is.na(ntrue), ntrue:=0]
  mer[,senscore:=ntrue/ncell]
  mer$cell=i
  myoutput <- rbind(myoutput, mer)
}
myoutput[,subject_id:=gsub("L|R", "", sample)]

## look at data
ggplot(myoutput, aes(senscore)) + geom_histogram() + facet_wrap(~cell, scales="free") 

## run linear models 
myoutput2 <- c()
for (i in names(scores)){
  cur <- myoutput[cell==i]
  cur <- cur[,list(senscore=mean(senscore)),subject_id] ##dealing with multiple samples per person
  img <- t(str_data_summarized[,cur$subject_id])
  for (j in colnames(img)){
    add <- unique(as.data.table( img[,j,drop=F], keep.rownames="subject_id"))
    colnames(add)[2] <- "imagingFeature"
    mer <- merge(cur, add)
    mer2 =as.data.frame(mer[,c("subject_id", "imagingFeature")])
    rownames(mer2) = mer$subject_id
    mer2$subject_id = NULL
    mer3 = t(mer2)
    cov = merge(mer, rownames(mer2)cov_str_summarized, by = "subject_id" )
    result <- as.data.table(summary(lm(senscore~imagingFeature, data=mer))$coef["imagingFeature",,drop=F])
    form = formula(paste("~senscore + (1|subject_id)"))
    voom_mer = voom(mer3)

    vobjDream <- voomWithDreamWeights(mer, form, cov, BPPARAM = param)

    dream(~senscore_percelltype + (1|subject_id), imagingfeature (voom object), metadata)
    dream(~senscore_percelltype + (1|IID) + age, imagingFeature, metadata)
    dream(~senscore_percelltype + (1|IID) + age + (1|sex), imagingFeature, metadata)

    result2 <- as.data.table(summary(lm(senscore~imagingFeature + IID, data=mer))$coef["imagingFeature",,drop=F])
    result3 <- as.data.table(summary(lm(senscore~imagingFeature + IID + age, data=mer))$coef["imagingFeature",,drop=F])
    result4 <- as.data.table(summary(lm(senscore~imagingFeature + IID + age + sex, data=mer))$coef["imagingFeature",,drop=F])

    result <- data.table(cell=i, imagingFeature=j, result)
    myoutput2 <- rbind(myoutput2,result)
  }
}
myoutput2[,padj:=p.adjust(`Pr(>|t|)`, "fdr")]


myoutput2[padj<0.05]
