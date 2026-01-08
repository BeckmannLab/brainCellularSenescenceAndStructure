### Load required libraries ----
library(batchtools)
library(devtools)
library(withr)
library(limma)
library(Glimma)
library(variancePartition)
library(sp)
library(biomaRt)
library(gsubfn)
library(data.table)
library(sp)
library(Matrix)
library(DropletUtils)
library(scater)
library(scRNAseq)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scDblFinder)
library(scran)
library(AnnotationHub)
library(BiocSingular)
library(BiocParallel)
library(uwot)
library(Seurat)
library(batchelor)
library(edgeR)
library(GSEABase)
library(AUCell)
library(bluster)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(muscat)
library(purrr)
library(UCell)
library(stringr)
library(rex)
library(data.table)
library(igraph)
library(foreach)
library(writexl)
library(doParallel)
registerDoParallel(cores=5)

##################
#set up
##################
### Load precomputed SenCor objects and define senescence gene sets ----
ls=load("/ssAnina_SenCor_Objects_wFishers_29JUL2024.RData") ##see scenic_fishers has the following three objects

# SenCorAll = the All nuclei Sen Score correlations both across all nuclei and separated by cell type. Also are fishers test SenFisher columns are the overlap of the 146 Sen Genes to each regulon target gene and the others are the POS or NEG DEGs for each cell type.
# SenCorMG = The MG only analysis
# SenCorExc= Exc only analysis

genes_42 <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", "CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", "GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", "LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", "MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", "STING1", "TGFB1", "TIMP2", "TNF", "TP53")
genes_125 <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")
SenGenes <- unique(c(genes_42, genes_125))

inAll=rownames(SenCorAll)

### Add Exc regulon annotations and export to Excel ----
SenCorExc$SenGene_Fisher_FDR = p.adjust(SenCorExc$SenGene_Fisher_p,method="fdr")
SenCorExc$geneName=gsub("(+)","",rownames(SenCorExc),fixed=TRUE)
intersect(SenGenes,SenCorExc$geneName)
SenCorExc$isSenGene = FALSE
SenCorExc$isSenGene[match(SenGenes,SenCorExc$geneName)[!is.na(match(SenGenes,SenCorExc$geneName))]]=TRUE
SenCorExc$inAll = FALSE
SenCorExc$inAll[match(inAll,rownames(SenCorExc))[!is.na(match(inAll,rownames(SenCorExc)))]]=TRUE
# SenCorExc[order(abs(SenCorExc$rho),decreasing=TRUE),grep("_p",colnames(SenCorExc),invert=TRUE)]
head(SenCorExc[order(abs(SenCorExc$rho),decreasing=TRUE),c("rho", "p", "Reg_nGenes", "SenGene_Fisher_OR", "SenGene_Fisher_FDR", "DEG_Fisher_OR", "DEG_Fisher_FDR", "POSDEG_Fisher_OR", "POSDEG_Fisher_FDR", "NEGDEG_Fisher_OR", "NEGDEG_Fisher_FDR","isSenGene","inAll")],20)
table(SenCorExc$inAll)

write_xlsx(SenCorExc, "/SenCorExc_regulons.xlsx")

### Add MG regulon annotations and export to Excel ----
SenCorMG$SenGene_Fisher_FDR = p.adjust(SenCorMG$SenGene_Fisher_p,method="fdr")
SenCorMG$geneName=gsub("(+)","",rownames(SenCorMG),fixed=TRUE)
intersect(SenGenes,SenCorMG$geneName)
SenCorMG$isSenGene = FALSE
SenCorMG$isSenGene[match(SenGenes,SenCorMG$geneName)[!is.na(match(SenGenes,SenCorMG$geneName))]]=TRUE
SenCorMG$inAll = FALSE
SenCorMG$inAll[match(inAll,rownames(SenCorMG))[!is.na(match(inAll,rownames(SenCorMG)))]]=TRUE
# SenCorMG[order(abs(SenCorMG$rho),decreasing=TRUE),grep("_p",colnames(SenCorMG),invert=TRUE)]
head(SenCorMG[order(abs(SenCorMG$rho),decreasing=TRUE),c("rho", "p", "Reg_nGenes", "SenGene_Fisher_OR", "SenGene_Fisher_FDR", "DEG_Fisher_OR", "DEG_Fisher_FDR", "POSDEG_Fisher_OR", "POSDEG_Fisher_FDR", "NEGDEG_Fisher_OR", "NEGDEG_Fisher_FDR","isSenGene","inAll")],20)
table(SenCorMG$inAll)

### Combine Exc + MG regulon metrics and export ----
SenCorExc$celltype = "Exc"
SenCorMG$celltype = "MG"

all = rbind(SenCorExc,SenCorMG )
all$regulon = rownames(all)

write_xlsx(all, "/SenCor_Exc_MG_regulons.xlsx")

### Reload core objects and re-run MG/Exc correlation checks ----
##load in 
excDE <- readRDS("/correct_edger_exc_5percent_fit.RDS")

mgDE <- readRDS("/correct_edger_mg_30percent_fit.RDS")

mg = readRDS("/MG_146_AUC_mat_top_30percent.RDS")
mg$cell_sen = "MG"

AUC <- as.data.frame(fread("/MG_auc_Final_ssAnina.csv"))

############
#MG
###########

##format
rownames(AUC) <- AUC$Regulon
AUC$Regulon <- NULL

#check
identical(rownames(mg), colnames(AUC))

AUCmg <- AUC

SenCorMG <- as.data.frame(matrix(nrow=nrow(AUCmg), ncol=2))
colnames(SenCorMG) <- c("rho", "p")
rownames(SenCorMG) <- rownames(AUCmg)

##check correlation
for (i in rownames(SenCorMG)){
  SenCorMG[i,]$rho <- cor.test(mg$sen_list, t(AUCmg)[,i],method="spearman")$estimate
        SenCorMG[i,]$p <- cor.test(mg$sen_list, t(AUCmg)[,i],method="spearman")$p.value

}

###########
#exc
##########

exc = readRDS("/Exc_146_AUC_mat_top_5percent.RDS")
exc$cell_sen = "Exc"

AUC <- as.data.frame(fread("/Exc_auc_Final_ssAnina.csv"))

# AUC <- getAUC(AUC)
rownames(AUC) <- AUC$Regulon
AUC$Regulon <- NULL

identical(rownames(exc), colnames(AUC))

AUCexc <- AUC

SenCorExc <- as.data.frame(matrix(nrow=nrow(AUCexc), ncol=2))
colnames(SenCorExc) <- c("rho", "p")
rownames(SenCorExc) <- rownames(AUCexc)

for (i in rownames(SenCorExc)){
  SenCorExc[i,]$rho <- cor.test(exc$sen_list, t(AUCexc)[,i],method="spearman")$estimate
        SenCorExc[i,]$p <- cor.test(exc$sen_list, t(AUCexc)[,i],method="spearman")$p.value

}

#####################################################################################
##########################
### Regulons and genes ###
############################

#####
##mg
######

reg <- fread(paste0("/MG_reg_Final_ssAnina.csv"))
x <- str_remove_all(reg$V9, rex(spaces))
all_matches <- str_match_all(x, rex("('", capture(not("'")), "',", capture(not(")")), ")"))
RegMG <- tstrsplit(rownames(AUCmg), split="(", fixed=T, keep=1)[[1]]

targetsMG <- list()
for (i in RegMG){
  index <- which(reg$V1==i)
  genes <- c()
   for (j in index){
   genes <- c(genes, all_matches[[j]][,2])
  }
 targetsMG[[i]] <- unique(genes)
} 

names(targetsMG) <- paste0(names(targetsMG),"(+)")

#######
###exc
#######
reg <- fread(paste0("/Exc_reg_Final_ssAnina.csv"))
x <- str_remove_all(reg$V9, rex(spaces))
all_matches <- str_match_all(x, rex("('", capture(not("'")), "',", capture(not(")")), ")"))
RegExc <- tstrsplit(rownames(AUCexc), split="(", fixed=T, keep=1)[[1]]


targetsExc <- list()
for (i in RegExc){
  index <- which(reg$V1==i)
  genes <- c()
   for (j in index){
   genes <- c(genes, all_matches[[j]][,2])
  }
 targetsExc[[i]] <- unique(genes)
} 

names(targetsExc) <- paste0(names(targetsExc),"(+)")
#####################################################################################

######################
# reload annotations
#########################
ls=load("/ssAnina_SenCor_Objects_wFishers_29JUL2024.RData")

# SenCorAll = the All nuclei Sen Score correlations both across all nuclei and separated by cell type. Also are fishers test SenFisher columns are the overlap of the 146 Sen Genes to each regulon target gene and the others are the POS or NEG DEGs for each cell type.
# SenCorMG = The MG only analysis
# SenCorExc= Exc only analysis

genes_42 <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", "CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", "GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", "LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", "MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", "STING1", "TGFB1", "TIMP2", "TNF", "TP53")
genes_125 <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")
SenGenes <- unique(c(genes_42, genes_125))

inAll=rownames(SenCorAll)

inAll2 = as.data.frame(inAll)

####
#exc
#####
SenCorExc$SenGene_Fisher_FDR = p.adjust(SenCorExc$SenGene_Fisher_p,method="fdr")
SenCorExc$rho_FDR = p.adjust(SenCorExc$p,method="fdr")
SenCorExc$geneName=gsub("(+)","",rownames(SenCorExc),fixed=TRUE)
intersect(SenGenes,SenCorExc$geneName)
SenCorExc$isSenGene = FALSE
SenCorExc$isSenGene[match(SenGenes,SenCorExc$geneName)[!is.na(match(SenGenes,SenCorExc$geneName))]]=TRUE
SenCorExc$inAll = FALSE
SenCorExc$inAll[match(inAll,rownames(SenCorExc))[!is.na(match(inAll,rownames(SenCorExc)))]]=TRUE
# SenCorExc[order(abs(SenCorExc$rho),decreasing=TRUE),grep("_p",colnames(SenCorExc),invert=TRUE)]
head(SenCorExc[order(abs(SenCorExc$rho),decreasing=TRUE),c("rho", "p", "Reg_nGenes", "SenGene_Fisher_OR", "SenGene_Fisher_FDR", "DEG_Fisher_OR", "DEG_Fisher_FDR", "POSDEG_Fisher_OR", "POSDEG_Fisher_FDR", "NEGDEG_Fisher_OR", "NEGDEG_Fisher_FDR","isSenGene","inAll")],20)
table(SenCorExc$inAll)


####
#MG
#####
SenCorMG$SenGene_Fisher_FDR = p.adjust(SenCorMG$SenGene_Fisher_p,method="fdr")
SenCorMG$rho_FDR = p.adjust(SenCorMG$p,method="fdr")
SenCorMG$geneName=gsub("(+)","",rownames(SenCorMG),fixed=TRUE)
intersect(SenGenes,SenCorMG$geneName)
SenCorMG$isSenGene = FALSE
SenCorMG$isSenGene[match(SenGenes,SenCorMG$geneName)[!is.na(match(SenGenes,SenCorMG$geneName))]]=TRUE
SenCorMG$inAll = FALSE
SenCorMG$inAll[match(inAll,rownames(SenCorMG))[!is.na(match(inAll,rownames(SenCorMG)))]]=TRUE
# SenCorMG[order(abs(SenCorMG$rho),decreasing=TRUE),grep("_p",colnames(SenCorMG),invert=TRUE)]
head(SenCorMG[order(abs(SenCorMG$rho),decreasing=TRUE),c("rho", "p", "Reg_nGenes", "SenGene_Fisher_OR", "SenGene_Fisher_FDR", "DEG_Fisher_OR", "DEG_Fisher_FDR", "POSDEG_Fisher_OR", "POSDEG_Fisher_FDR", "NEGDEG_Fisher_OR", "NEGDEG_Fisher_FDR","isSenGene","inAll")],20)
table(SenCorMG$inAll)

#####################################################################################

################
# create network
################

#########
#mg
##########
### Build MG regulon–regulon overlap network (graph edges) ----
regulons=targetsMG
allGenes=rownames(mgDE)
names(regulons)=gsub("(+)","",names(regulons),fixed=TRUE)

parents=list()
for(i in 1:length(names(regulons))){
  regName=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName %in% x }))))
  parentsVec=setdiff(parentsVec,regName)
  if(length(parentsVec)>0){
    parents[[regName]]=data.frame(child=regName, parent=parentsVec)
  }
}
parents_df=do.call(rbind,parents)


parents=list()
for(i in 1:length(names(regulons))){
  regName1=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName1 %in% x }))))
  parentsVec=setdiff(parentsVec,regName1)
  if(length(parentsVec)>0){
    parents[[regName1]]=data.frame(parent=parentsVec, child=regName1, p.val=NA, OR=NA, inReg1_inReg2=NA, notInReg1_inReg2=NA, inReg1_notInReg2=NA, notInReg1_notInReg2=NA)
    for(j in 1:length(parentsVec)){
      regName2=parentsVec[j]
      allGenesInRegs=data.frame(geneName=allGenes)
      allGenesInRegs$inReg1="notInReg1"
      allGenesInRegs$inReg1[allGenesInRegs$geneName %in% regulons[[regName1]]]="inReg1"
      allGenesInRegs$inReg1=factor(allGenesInRegs$inReg1,levels=c("inReg1", "notInReg1"))
      allGenesInRegs$inReg2="notInReg2"
      allGenesInRegs$inReg2[allGenesInRegs$geneName %in% regulons[[regName2]]]="inReg2"
      allGenesInRegs$inReg2=factor(allGenesInRegs$inReg2,levels=c("inReg2", "notInReg2"))
      current=table(allGenesInRegs$inReg1,allGenesInRegs$inReg2)
      test=fisher.test(current,alternative="greater")
      parents[[regName1]]$p.val[j]=test$p.value
      parents[[regName1]]$OR[j]=test$estimate
      parents[[regName1]]$inReg1_inReg2[j]=current["inReg1","inReg2"]
      parents[[regName1]]$notInReg1_inReg2[j]=current["notInReg1","inReg2"]
      parents[[regName1]]$inReg1_notInReg2[j]=current["inReg1","notInReg2"]
      parents[[regName1]]$notInReg1_notInReg2[j]=current["notInReg1","notInReg2"]

    }
  }
}
parents_df=do.call(rbind,parents)
parents_df$adj.p.val=p.adjust(parents_df$p.val, method="bonferroni")
parents_df[order(parents_df$inReg1_inReg2,decreasing=TRUE),]
realEdges=parents_df[parents_df$adj.p.val<0.05,]
to_keep=SenCorMG$geneName[SenCorMG$rho_FDR<0.05 & (SenCorMG$POSDEG_Fisher_FDR<0.05 |SenCorMG$NEGDEG_Fisher_FDR<0.05)]
to_keep=intersect(intersect(realEdges$parent,realEdges$child),to_keep)
identical(sort(unique(realEdges$parent[realEdges$parent %in% to_keep])),sort(to_keep))
identical(sort(unique(realEdges$child[realEdges$child %in% to_keep])),sort(to_keep))
realEdges=realEdges[realEdges$parent %in% to_keep & realEdges$child %in% to_keep,]
# intersect(intersect(realEdges$parent,realEdges$child),to_keep)
intersect(realEdges$parent,realEdges$child)
# fwrite(realEdges,"~/www/results/mg_graph.txt")


###########
#exc
###########

### Build Exc regulon–regulon overlap skeleton (graph edges) ----
regulons=targetsExc
allGenes=rownames(excDE)
names(regulons)=gsub("(+)","",names(regulons),fixed=TRUE)

parents=list()
for(i in 1:length(names(regulons))){
  regName=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName %in% x }))))
  parentsVec=setdiff(parentsVec,regName)
  if(length(parentsVec)>0){
    parents[[regName]]=data.frame(child=regName, parent=parentsVec)
  }
}
parents_df=do.call(rbind,parents)

################################################################
##see which regulon enriched for imaging DE and subset to those 
##################################################################

#load in imaging DE list
myres = readRDS("/single_DE_results_6.13.24.RDS")
myres$feature = gsub("scale", "", myres$feature)
rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")
# Remove rows with specified row names
myres2 <- myres[!myres$feature %in% rows_to_remove, ]

myres3 = myres2[which(myres2$assay == "MG"),]

features = unique(myres3$feature)
cell_type = unique(myres3$assay)

resFinal = c()
foreach (x = features) %do% {
    df = myres3[which(myres3$feature == x),]
    df=df[order(df$P.Value, decreasing = FALSE),]
    df$DE = "notSignif"
    cutoff = floor((1-limma::propTrueNull(df$P.Value))*nrow(df))
    df$DE[1:cutoff] = "Signif"
    resFinal = rbind(df, resFinal)
}

####
#mg
####
### MG: test regulons for enrichment of imaging DEGs (sMRI features) ----
regulons=targetsMG
allGenes=rownames(mgDE)
names(regulons)=gsub("(+)","",names(regulons),fixed=TRUE)

parents=list()
for(i in 1:length(names(regulons))){
  regName=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName %in% x }))))
  parentsVec=setdiff(parentsVec,regName)
  if(length(parentsVec)>0){
    parents[[regName]]=data.frame(child=regName, parent=parentsVec)
  }
}
parents_df=do.call(rbind,parents)

resFinal$de=paste0(resFinal$feature,resFinal$assay)

parents=c()

res=foreach(i = 1:length(names(regulons))) %do% { #for each regulon
  print(i)
  regName1=names(regulons)[i] ##select regulon i parent
  # parentsVec=names(which(unlist(lapply(regulons, function(x){ regName1 %in% x })))) #looking to see if reName1 xist in other regulons
  # parentsVec=setdiff(parentsVec,regName1) #removing regulon controlled by that gene
  inRegulon = regulons[[regName1]]
  if(length(inRegulon)>0){ #if there are genes in the regulon 
    grid=data.frame(expand.grid(unique(resFinal$feature),unique(resFinal$assay)))
    colnames(grid)=c("sMRI_feature","cellType")
    parents=cbind(grid,data.frame(parent=regName1, p.val=NA, OR=NA, inReg1_DE=NA, notInReg1_DE=NA, inReg1_notDE=NA, notInReg1_notDE=NA))
    for(j in 1:nrow(grid)){ #for each gene in the regulon (i)
      subset=resFinal[resFinal$de==paste0(grid$sMRI_feature[j],grid$cellType[j]),]
      smriName=unique(subset$feature) #imaging
      cellType=unique(subset$assay)
      allGenesInRegs=data.frame(geneName=allGenes) #all genes in regulon
      allGenesInRegs$inReg1="notInReg1" #set up to define which ones are in reg1
      allGenesInRegs$inReg1[allGenesInRegs$geneName %in% inRegulon]="inReg1" ##which ones are in regulon 1
      allGenesInRegs$inReg1=factor(allGenesInRegs$inReg1,levels=c("inReg1", "notInReg1")) #factor and set levels
      allGenesInRegs$deSmri="notDESmri" #not in smri 
      allGenesInRegs$deSmri[allGenesInRegs$geneName %in% subset$ID[subset$DE=="Signif"]]="deSmri" #which ones are in smri
      allGenesInRegs$deSmri=factor(allGenesInRegs$deSmri,levels=c("deSmri", "notDESmri")) #factor and set levels
      current=table(allGenesInRegs$inReg1,allGenesInRegs$deSmri) #table between in smri and regulon 1
      test=fisher.test(current,alternative="greater") #fishers test
      print(j)
      parents$p.val[j]=test$p.value #define p value
      parents$OR[j]=test$estimate #rho 
      parents$inReg1_DE[j]=current["inReg1","deSmri"] #rename 
      parents$notInReg1_DE[j]=current["notInReg1","deSmri"] #rename 
      parents$inReg1_notDE[j]=current["inReg1","notDESmri"] #rename 
      parents$notInReg1_notDE[j]=current["notInReg1","notDESmri"] #rename 
      parents$sMRI_feature[j]=smriName
      parents$cellType[j]=cellType

    }
  }
  parents
}


res_df=do.call(rbind,res)
res_df$adj.p.val=p.adjust(res_df$p.val, method="bonferroni")

nodes_withSMRIDE=unique(res_df$parent[res_df$adj.p.val<0.05])


### Rebuild MG regulon–regulon network and filter by sMRI-enriched regulons ----
regulons=targetsMG
allGenes=rownames(mgDE)
names(regulons)=gsub("(+)","",names(regulons),fixed=TRUE)

parents=list()
for(i in 1:length(names(regulons))){
  regName=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName %in% x }))))
  parentsVec=setdiff(parentsVec,regName)
  if(length(parentsVec)>0){
    parents[[regName]]=data.frame(child=regName, parent=parentsVec)
  }
}
parents_df=do.call(rbind,parents)


parents=list()
for(i in 1:length(names(regulons))){
  regName1=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName1 %in% x }))))
  parentsVec=setdiff(parentsVec,regName1)
  if(length(parentsVec)>0){
    parents[[regName1]]=data.frame(parent=parentsVec, child=regName1, p.val=NA, OR=NA, inReg1_inReg2=NA, notInReg1_inReg2=NA, inReg1_notInReg2=NA, notInReg1_notInReg2=NA)
    for(j in 1:length(parentsVec)){
      regName2=parentsVec[j]
      allGenesInRegs=data.frame(geneName=allGenes)
      allGenesInRegs$inReg1="notInReg1"
      allGenesInRegs$inReg1[allGenesInRegs$geneName %in% regulons[[regName1]]]="inReg1"
      allGenesInRegs$inReg1=factor(allGenesInRegs$inReg1,levels=c("inReg1", "notInReg1"))
      allGenesInRegs$inReg2="notInReg2"
      allGenesInRegs$inReg2[allGenesInRegs$geneName %in% regulons[[regName2]]]="inReg2"
      allGenesInRegs$inReg2=factor(allGenesInRegs$inReg2,levels=c("inReg2", "notInReg2"))
      current=table(allGenesInRegs$inReg1,allGenesInRegs$inReg2)
      test=fisher.test(current,alternative="greater")
      parents[[regName1]]$p.val[j]=test$p.value
      parents[[regName1]]$OR[j]=test$estimate
      parents[[regName1]]$inReg1_inReg2[j]=current["inReg1","inReg2"]
      parents[[regName1]]$notInReg1_inReg2[j]=current["notInReg1","inReg2"]
      parents[[regName1]]$inReg1_notInReg2[j]=current["inReg1","notInReg2"]
      parents[[regName1]]$notInReg1_notInReg2[j]=current["notInReg1","notInReg2"]

    }
  }
}
parents_df=do.call(rbind,parents)
parents_df$adj.p.val=p.adjust(parents_df$p.val, method="bonferroni")
parents_df[order(parents_df$inReg1_inReg2,decreasing=TRUE),]
realEdges=parents_df[parents_df$adj.p.val<0.05,]
to_keep=SenCorMG$geneName[SenCorMG$rho_FDR<0.05 & (SenCorMG$POSDEG_Fisher_FDR<0.05 |SenCorMG$NEGDEG_Fisher_FDR<0.05)]
to_keep=intersect(to_keep, nodes_withSMRIDE) ##only keep ones that have DE for sMRI

to_keep=intersect(intersect(realEdges$parent,realEdges$child),to_keep)
identical(sort(unique(realEdges$parent[realEdges$parent %in% to_keep])),sort(to_keep))
identical(sort(unique(realEdges$child[realEdges$child %in% to_keep])),sort(to_keep))
realEdges=realEdges[realEdges$parent %in% to_keep & realEdges$child %in% to_keep,]
# intersect(intersect(realEdges$parent,realEdges$child),to_keep)
intersect(realEdges$parent,realEdges$child)
#  [1] "ETV6"   "FLI1"   "JAZF1"  "RUNX1"  "ELF1"   "IKZF1"  "ZNF148" "FOXN2" 
#  [9] "HBP1"   "IRF8"   "SPI1"   "FOXP2"  "ELK3"   "ETV5"   "FOXN3"  "CREB5" 
# [17] "MXI1"   "CREB1"  "REL"    "TAL1"   "ARID3A"

###########
#exc
###########

### Exc
#load in imaging DE list
myres = readRDS("single_DE_results_6.13.24.RDS")
myres$feature = gsub("scale", "", myres$feature)
rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")
# Remove rows with specified row names
myres2 <- myres[!myres$feature %in% rows_to_remove, ]

myres3 = myres2[which(myres2$assay == "MG"),]

features = unique(myres3$feature)
cell_type = unique(myres3$assay)

resFinal = c()
foreach (x = features) %do% {
    df = myres3[which(myres3$feature == x),]
    df=df[order(df$P.Value, decreasing = FALSE),]
    df$DE = "notSignif"
    cutoff = floor((1-limma::propTrueNull(df$P.Value))*nrow(df))
    df$DE[1:cutoff] = "Signif"
    resFinal = rbind(df, resFinal)
   }


regulons=targetsExc
allGenes=rownames(excDE)
names(regulons)=gsub("(+)","",names(regulons),fixed=TRUE)

parents=list()
for(i in 1:length(names(regulons))){
  regName=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName %in% x }))))
  parentsVec=setdiff(parentsVec,regName)
  if(length(parentsVec)>0){
    parents[[regName]]=data.frame(child=regName, parent=parentsVec)
  }
}
parents_df=do.call(rbind,parents)

resFinal$de=paste0(resFinal$feature,resFinal$assay)

parents=c()


res=foreach(i = 1:length(names(regulons))) %dopar% { #for each regulon
  print(i)
  regName1=names(regulons)[i] ##select regulon i parent
  # parentsVec=names(which(unlist(lapply(regulons, function(x){ regName1 %in% x })))) #looking to see if reName1 xist in other regulons
  # parentsVec=setdiff(parentsVec,regName1) #removing regulon controlled by that gene
  inRegulon = regulons[[regName1]]
  if(length(inRegulon)>0){ #if there are genes in the regulon 
    grid=data.frame(expand.grid(unique(resFinal$feature),unique(resFinal$assay)))
    colnames(grid)=c("sMRI_feature","cellType")
    parents=cbind(grid,data.frame(parent=regName1, p.val=NA, OR=NA, inReg1_DE=NA, notInReg1_DE=NA, inReg1_notDE=NA, notInReg1_notDE=NA))
    for(j in 1:nrow(grid)){ #for each gene in the regulon (i)
      subset=resFinal[resFinal$de==paste0(grid$sMRI_feature[j],grid$cellType[j]),]
      smriName=unique(subset$feature) #imaging
      cellType=unique(subset$assay)
      allGenesInRegs=data.frame(geneName=allGenes) #all genes in regulon
      allGenesInRegs$inReg1="notInReg1" #set up to define which ones are in reg1
      allGenesInRegs$inReg1[allGenesInRegs$geneName %in% inRegulon]="inReg1" ##which ones are in regulon 1
      allGenesInRegs$inReg1=factor(allGenesInRegs$inReg1,levels=c("inReg1", "notInReg1")) #factor and set levels
      allGenesInRegs$deSmri="notDESmri" #not in smri 
      allGenesInRegs$deSmri[allGenesInRegs$geneName %in% subset$ID[subset$DE=="Signif"]]="deSmri" #which ones are in smri
      allGenesInRegs$deSmri=factor(allGenesInRegs$deSmri,levels=c("deSmri", "notDESmri")) #factor and set levels
      current=table(allGenesInRegs$inReg1,allGenesInRegs$deSmri) #table between in smri and regulon 1
      test=fisher.test(current,alternative="greater") #fishers test
      print(j)
      parents$p.val[j]=test$p.value #define p value
      parents$OR[j]=test$estimate #rho 
      parents$inReg1_DE[j]=current["inReg1","deSmri"] #rename 
      parents$notInReg1_DE[j]=current["notInReg1","deSmri"] #rename 
      parents$inReg1_notDE[j]=current["inReg1","notDESmri"] #rename 
      parents$notInReg1_notDE[j]=current["notInReg1","notDESmri"] #rename 
      parents$sMRI_feature[j]=smriName
      parents$cellType[j]=cellType

    }
  }
  parents
}

res_df=do.call(rbind,res)
res_df$adj.p.val=p.adjust(res_df$p.val, method="bonferroni")

nodes_withSMRIDE=unique(res_df$parent[res_df$adj.p.val<0.05])

# [1] "BCLAF1"  "BHLHE40" "BPTF"    "BRF1"    "CBFB"    "CUX1"    "CUX2"   
#  [8] "ELF1"    "ELF2"    "FOXN2"   "FOXN3"   "FOXO1"   "FOXP1"   "FOXP2"  
# [15] "HSF4"    "JUN"     "JUND"    "NFIC"    "NFIX"    "PBX1"    "PKNOX2" 
# [22] "POU3F3"  "PPARA"   "RAD21"   "RBPJ"    "RXRA"    "SMARCA4" "SOX10"  
# [29] "SOX2"    "SOX8"    "SREBF1"  "SREBF2"  "VEZF1"   "ZEB1"    "ZIC1"   
# [36] "ZMIZ1"   "ZNF226"  "ZNF91"  

### Build Exc regulon-regulon network and filter by sMRI-enriched regulons ----
parents=list()
for(i in 1:length(names(regulons))){
  regName1=names(regulons)[i]
  parentsVec=names(which(unlist(lapply(regulons, function(x){ regName1 %in% x }))))
  parentsVec=setdiff(parentsVec,regName1)
  if(length(parentsVec)>0){
    parents[[regName1]]=data.frame(parent=parentsVec, child=regName1, p.val=NA, OR=NA, inReg1_inReg2=NA, notInReg1_inReg2=NA, inReg1_notInReg2=NA, notInReg1_notInReg2=NA)
    for(j in 1:length(parentsVec)){
      regName2=parentsVec[j]
      allGenesInRegs=data.frame(geneName=allGenes)
      allGenesInRegs$inReg1="notInReg1"
      allGenesInRegs$inReg1[allGenesInRegs$geneName %in% regulons[[regName1]]]="inReg1"
      allGenesInRegs$inReg1=factor(allGenesInRegs$inReg1,levels=c("inReg1", "notInReg1"))
      allGenesInRegs$inReg2="notInReg2"
      allGenesInRegs$inReg2[allGenesInRegs$geneName %in% regulons[[regName2]]]="inReg2"
      allGenesInRegs$inReg2=factor(allGenesInRegs$inReg2,levels=c("inReg2", "notInReg2"))
      current=table(allGenesInRegs$inReg1,allGenesInRegs$inReg2)
      test=fisher.test(current,alternative="greater")
      parents[[regName1]]$p.val[j]=test$p.value
      parents[[regName1]]$OR[j]=test$estimate
      parents[[regName1]]$inReg1_inReg2[j]=current["inReg1","inReg2"]
      parents[[regName1]]$notInReg1_inReg2[j]=current["notInReg1","inReg2"]
      parents[[regName1]]$inReg1_notInReg2[j]=current["inReg1","notInReg2"]
      parents[[regName1]]$notInReg1_notInReg2[j]=current["notInReg1","notInReg2"]

    }
  }
}

parents_df=do.call(rbind,parents)
parents_df$adj.p.val=p.adjust(parents_df$p.val, method="fdr")
parents_df[order(parents_df$inReg1_inReg2,decreasing=TRUE),]
realEdges=parents_df[parents_df$adj.p.val<0.05,]
to_keep=SenCorExc$geneName[SenCorExc$rho_FDR<0.05 & (SenCorExc$POSDEG_Fisher_FDR<0.05 |SenCorExc$NEGDEG_Fisher_FDR<0.05)]
to_keep=intersect(to_keep, nodes_withSMRIDE) ##only keep ones that have DE for sMRI


to_keep=intersect(intersect(realEdges$parent,realEdges$child),to_keep)
identical(sort(unique(realEdges$parent[realEdges$parent %in% to_keep])),sort(to_keep))
identical(sort(unique(realEdges$child[realEdges$child %in% to_keep])),sort(to_keep))
realEdges=realEdges[realEdges$parent %in% to_keep & realEdges$child %in% to_keep,]

intersect(realEdges$parent,realEdges$child)
