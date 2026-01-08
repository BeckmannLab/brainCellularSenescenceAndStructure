### Load required packages ----
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

### Microglia (MG): load AUC matrix and correlate regulon activity with senescence score ----
mg = readRDS("/MG_146_AUC_mat_top_30percent.RDS") #microglia 
mg$cell_sen = "MG"

AUC <- as.data.frame(fread("/MG_auc_Final_ssAnina.csv"))
rownames(AUC) <- AUC$Regulon
AUC$Regulon <- NULL

identical(rownames(mg), colnames(AUC))

AUCmg <- AUC

SenCorMG <- as.data.frame(matrix(nrow=nrow(AUCmg), ncol=2))
colnames(SenCorMG) <- c("rho", "p")
rownames(SenCorMG) <- rownames(AUCmg)

for (i in rownames(SenCorMG)){
	SenCorMG[i,]$rho <- cor.test(mg$sen_list, t(AUCmg)[,i],method="spearman")$estimate
        SenCorMG[i,]$p <- cor.test(mg$sen_list, t(AUCmg)[,i],method="spearman")$p.value

}

### Excitatory neurons (Exc): load AUC matrix and correlate regulon activity with senescence score ----
exc = readRDS("/Exc_146_AUC_mat_top_5percent.RDS")
exc$cell_sen = "Exc"

AUC <- as.data.frame(fread("/Exc_auc_Final_ssAnina.csv"))

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


### Regulons and genes ###

### Derive target genes for MG regulons ----
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


### Derive target genes for Exc regulons ----
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

### Define senescence-associated gene sets (42 and 125 markers) ----
genes_42 <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", 
"CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", 
"GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", 
"LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", 
"MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", 
"STING1", "TGFB1", "TIMP2", "TNF", "TP53")


genes_125 <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
"BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", 
"CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", 
"CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", 
"CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", 
"CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", 
"FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", 
"HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", 
"IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", 
"IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", 
"INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", 
"MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", 
"NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", 
"PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", 
"SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", 
"TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", 
"TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")

SenGenes <- unique(c(genes_42, genes_125))

### Load cell-type-specific DE results for Exc and MG and annotate with senescence / regulon targets ----
excDE <- readRDS("/correct_edger_exc_5percent_fit.RDS")

mgDE <- readRDS("/correct_edger_mg_30percent_fit.RDS")


excDE$SenGenes <- rownames(excDE) %in% SenGenes

for (i in names(targetsExc)){
        excDE$DEG <- excDE$FDR<=0.05
        excDE$POSDEG <- excDE$FDR<=0.05 & excDE$logFC>0
	excDE$NEGDEG <- excDE$FDR<=0.05 & excDE$logFC<=0
        excDE[,i] <- rownames(excDE) %in% targetsExc[[i]]
}

mgDE$SenGenes <- rownames(mgDE) %in% SenGenes

for (i in names(targetsMG)){
	mgDE$DEG <- mgDE$FDR<=0.05
        mgDE$POSDEG <- mgDE$FDR<=0.05 & mgDE$logFC>0
	mgDE$NEGDEG <- mgDE$FDR<=0.05 & mgDE$logFC<=0
	mgDE[,i] <- rownames(mgDE) %in% targetsMG[[i]]
}

### MG: Fisher enrichment tests for regulon targets vs senescence / DE gene sets ----
SenCorMG$Reg_nGenes <- NA
SenCorMG$SenGene_Fisher_OR <- NA
SenCorMG$SenGene_Fisher_p <- NA 
SenCorMG$DEG_Fisher_OR <- NA
SenCorMG$DEG_Fisher_p <- NA
SenCorMG$POSDEG_Fisher_OR <- NA 
SenCorMG$POSDEG_Fisher_p <- NA 
SenCorMG$NEGDEG_Fisher_OR <- NA 
SenCorMG$NEGDEG_Fisher_p <- NA 


for (i in names(targetsMG)){
	SenCorMG[i,]$Reg_nGenes <- length(targetsMG[[i]])
	SenCorMG[i,]$SenGene_Fisher_OR <- fisher.test(mgDE$SenGenes, mgDE[,i])$estimate
	SenCorMG[i,]$SenGene_Fisher_p <- fisher.test(mgDE$SenGenes, mgDE[,i])$p.value

	SenCorMG[i,]$DEG_Fisher_OR <- fisher.test(mgDE$DEG, mgDE[,i])$estimate
	SenCorMG[i,]$DEG_Fisher_p <- fisher.test(mgDE$DEG, mgDE[,i])$p.value

	SenCorMG[i,]$POSDEG_Fisher_OR <- fisher.test(mgDE$POSDEG, mgDE[,i])$estimate
	SenCorMG[i,]$POSDEG_Fisher_p <- fisher.test(mgDE$POSDEG, mgDE[,i])$p.value

	SenCorMG[i,]$NEGDEG_Fisher_OR <- fisher.test(mgDE$NEGDEG, mgDE[,i])$estimate
	SenCorMG[i,]$NEGDEG_Fisher_p <- fisher.test(mgDE$NEGDEG, mgDE[,i])$p.value

}

SenCorMG$DEG_Fisher_FDR <- p.adjust(SenCorMG$DEG_Fisher_p, method="fdr")
SenCorMG$POSDEG_Fisher_FDR <- p.adjust(SenCorMG$POSDEG_Fisher_p, method="fdr")
SenCorMG$NEGDEG_Fisher_FDR <- p.adjust(SenCorMG$NEGDEG_Fisher_p, method="fdr")

### Exc: Fisher enrichment tests for regulon targets vs senescence / DE gene sets ----
SenCorExc$Reg_nGenes <- NA
SenCorExc$SenGene_Fisher_OR <- NA
SenCorExc$SenGene_Fisher_p <- NA 
SenCorExc$DEG_Fisher_OR <- NA
SenCorExc$DEG_Fisher_p <- NA
SenCorExc$POSDEG_Fisher_OR <- NA 
SenCorExc$POSDEG_Fisher_p <- NA 
SenCorExc$NEGDEG_Fisher_OR <- NA 
SenCorExc$NEGDEG_Fisher_p <- NA 

for (i in names(targetsExc)){
	SenCorExc[i,]$Reg_nGenes <- length(targetsExc[[i]])
	SenCorExc[i,]$SenGene_Fisher_OR <- fisher.test(excDE$SenGenes, excDE[,i])$estimate
	SenCorExc[i,]$SenGene_Fisher_p <- fisher.test(excDE$SenGenes, excDE[,i])$p.value

	SenCorExc[i,]$DEG_Fisher_OR <- fisher.test(excDE[,i], excDE$DEG)$estimate
	SenCorExc[i,]$DEG_Fisher_p <- fisher.test(excDE[,i], excDE$DEG)$p.value

	SenCorExc[i,]$POSDEG_Fisher_OR <- fisher.test(excDE[,i], excDE$POSDEG)$estimate
	SenCorExc[i,]$POSDEG_Fisher_p <- fisher.test(excDE[,i], excDE$POSDEG)$p.value

	SenCorExc[i,]$NEGDEG_Fisher_OR <- fisher.test(excDE[,i], excDE$NEGDEG)$estimate
	SenCorExc[i,]$NEGDEG_Fisher_p <- fisher.test(excDE[,i], excDE$NEGDEG)$p.value

}

SenCorExc$DEG_Fisher_FDR <- p.adjust(SenCorExc$DEG_Fisher_p, method="fdr")
SenCorExc$POSDEG_Fisher_FDR <- p.adjust(SenCorExc$POSDEG_Fisher_p, method="fdr")
SenCorExc$NEGDEG_Fisher_FDR <- p.adjust(SenCorExc$NEGDEG_Fisher_p, method="fdr")

### All Nuclei Analusis ###
### Compile senescence scores across all major cell types and merge with global AUC matrix ----
sen <- list()

sen$mg = readRDS("/MG_146_AUC_mat_top_30percent.RDS")
sen$mg$cell_sen = "MG"

sen$exc = readRDS("/Exc_146_AUC_mat_top_5percent.RDS")
sen$exc$cell_sen = "Exc"

sen$oli = readRDS("/Oli_146_AUC_mat_top_10percent.RDS")
sen$oli$cell_sen = "Oli"

sen$int = readRDS("/Int_146_AUC_mat_top_5percent.RDS")
sen$int$cell_sen = "Int"

sen$noneu = readRDS("/NonNeu_146_AUC_mat_top_1percent.RDS")
sen$noneu$cell_sen = "NonNeu"

sen$ast = readRDS("/Ast_146_AUC_mat_top_1percent.RDS")
sen$ast$cell_sen = "Ast"

sen$opc = readRDS("/OPC_146_AUC_mat_top_1percent.RDS")
sen$opc$cell_sen = "OPC"

for (i in names(sen)){
 sen[[i]]$cell <- rownames(sen[[i]])
}

senAll <- do.call(rbind, sen)

AUC <- as.data.frame(fread("/All_auc_Final_ssAnina.csv"))  ##output from SCENIC

rownames(AUC) <- AUC$Regulon
AUC$Regulon <- NULL

AUCall <- as.data.frame(t(AUC))
AUCall$cell <- rownames(AUCall)

senMerge <- merge(senAll, AUCall, by="cell")

### Compute regulon–senescence correlations per cell class and across all nuclei ----
SenCorAll <- as.data.frame(matrix(nrow=nrow(AUC), ncol=length(unique(senMerge$cell_sen))))
colnames(SenCorAll) <- unique(senMerge$cell_sen)
rownames(SenCorAll) <- rownames(AUC)

for (i in colnames(SenCorAll)){
 for (j in rownames(SenCorAll)){
	SenCorAll[j,i] <- cor.test(senMerge[senMerge$cell_sen==i,]$sen_list, senMerge[senMerge$cell_sen==i,j],method="spearman")$estimate
 }
}


### Fisher Test ###
# SenGene Overlap #

### Build regulon target lists for All-nuclei SCENIC run and run gene-level senescence Fisher tests ----
load("/ExpressedGenes_for_ssAnina_Scenic.RData") #expressed genes

reg <- fread(paste0("/All_reg_Final_ssAnina.csv"))
x <- str_remove_all(reg$V9, rex(spaces))
all_matches <- str_match_all(x, rex("('", capture(not("'")), "',", capture(not(")")), ")"))
Reg <- tstrsplit(rownames(AUC), split="(", fixed=T, keep=1)[[1]]


targetsAll <- list()
for (i in Reg){
  index <- which(reg$V1==i)
  genes <- c()
   for (j in index){
   genes <- c(genes, all_matches[[j]][,2])
  }
 targetsAll[[i]] <- unique(genes)
} 

names(targetsAll) <- paste0(names(targetsAll),"(+)")


SenCorAll$AllNuc <- NA
SenCorAll$nGenes <- NA

for (i in rownames(SenCorAll)){
	SenCorAll[i,]$nGenes <- length(targetsAll[[i]])
	SenCorAll[i,]$AllNuc <- cor.test(senMerge$sen_list, senMerge[,i],method="spearman")$estimate
}


SenFisher <- as.data.frame(matrix(ncol=1, nrow=length(genes_all)))
rownames(SenFisher) <- genes_all
SenFisher$sen_list <- rownames(SenFisher) %in% SenGenes
for (i in names(targetsAll)){
	SenFisher[,i] <- rownames(SenFisher) %in% targetsAll[[i]]
}

SenCorAll$SenFisherOR <- NA
SenCorAll$SenFisherP <- NA
for (i in names(targetsAll)){
	SenCorAll[i,]$SenFisherOR <- fisher.test(SenFisher$sen_list, SenFisher[,i])$estimate
	SenCorAll[i,]$SenFisherP <- fisher.test(SenFisher$sen_list, SenFisher[,i])$p.value
}

SenCorAll$SenFisherFDR <- p.adjust(SenCorAll$SenFisherP, method="fdr")


# SenDE overlap #
### Load DE results for all cell types and annotate DEG status ----
ssDE <- list()

ssDE$Exc = readRDS("/correct_edger_exc_5percent_fit.RDS")

ssDE$Oli = readRDS("/correct_edger_oli_10percent_fit.RDS")

ssDE$Int = readRDS("/correct_edger_int_5percent_fit.RDS")

ssDE$NonNeu = readRDS("/correct_edger_noneu_1percent_fit.RDS")

ssDE$MG = readRDS("/correct_edger_mg_30percent_fit.RDS")

ssDE$Ast = readRDS("/correct_edger_ast_1percent_fit.RDS")

ssDE$OPC = readRDS("/correct_edger_opc_1percent_fit.RDS")

for (i in names(ssDE)){
        ssDE[[i]]$DEG <- ssDE[[i]]$FDR<=0.05
        ssDE[[i]]$POSDEG <- ssDE[[i]]$FDR<=0.05 & ssDE[[i]]$logFC>0
	ssDE[[i]]$NEGDEG <- ssDE[[i]]$FDR<=0.05 & ssDE[[i]]$logFC<=0
}

### Build DEG membership matrix (per cell type, POS/NEG) for senescence Fisher tests ----
for (i in names(ssDE)){
	SenFisher[,paste0(i,"_POSDEG")] <- rownames(SenFisher) %in% rownames(ssDE[[i]][ssDE[[i]]$POSDEG=="TRUE",])
	SenFisher[,paste0(i,"_NEGDEG")] <- rownames(SenFisher) %in% rownames(ssDE[[i]][ssDE[[i]]$NEGDEG=="TRUE",])
}

### Regulon–DEG enrichment Fisher tests across all cell types ----
for (i in names(ssDE)){
	SenCorAll[,paste0(i,"_POSDEG_FisherOR")] <- NA
	SenCorAll[,paste0(i,"_POSDEG_FisherP")] <- NA

	SenCorAll[,paste0(i,"_NEGDEG_FisherOR")] <- NA
	SenCorAll[,paste0(i,"_NEGDEG_FisherP")] <- NA

	for (j in names(targetsAll)){
		if(length(which(SenFisher[,paste0(i,"_POSDEG")]==TRUE)) > 0){
		 SenCorAll[j,paste0(i,"_POSDEG_FisherOR")] <- fisher.test(SenFisher[,j], SenFisher[,paste0(i,"_POSDEG")])$estimate
		 SenCorAll[j,paste0(i,"_POSDEG_FisherP")] <- fisher.test(SenFisher[,j], SenFisher[,paste0(i,"_POSDEG")])$p.value
		}

		if(length(which(SenFisher[,paste0(i,"_NEGDEG")]==TRUE)) > 0){
		 SenCorAll[j,paste0(i,"_NEGDEG_FisherOR")] <- fisher.test(SenFisher[,j], SenFisher[,paste0(i,"_NEGDEG")])$estimate
		 SenCorAll[j,paste0(i,"_NEGDEG_FisherP")] <- fisher.test(SenFisher[,j], SenFisher[,paste0(i,"_NEGDEG")])$p.value
		}
	} 
	SenCorAll[,paste0(i,"_POSDEG_FisherFDR")] <- p.adjust(SenCorAll[,paste0(i,"_POSDEG_FisherP")], method="fdr")
	if(length(which(SenFisher[,paste0(i,"_NEGDEG")]==TRUE)) > 0){
	 SenCorAll[,paste0(i,"_NEGDEG_FisherFDR")] <- p.adjust(SenCorAll[,paste0(i,"_NEGDEG_FisherP")], method="fdr")
	}
	
}

### Save all correlation and Fisher test results ----
save(SenCorAll, SenCorMG, SenCorExc, file = "ssAnina_SenCor_Objects_wFishers_29JUL2024.RData")
