library(qs)
library(SingleCellExperiment)
library(Seurat)
library(foreach)
library(doParallel)
registerDoParallel(cores=40)


cmc = qread("sce_all_CMC.qs") #see setup

#subset by cell type 
celltypes = c("Exc", "Int", "Astro", "Oligo", "OPC", "Micro")
resfinal=foreach(x = celltypes) %dopar% {
	print(x)
	sub_cmc = cmc[, which(colData(cmc)$new_celltype == x)]
	df_cmc =counts(sub_cmc)
	df_cmc = as.matrix(df_cmc)
	df_cmc_t = t(df_cmc)
	write.csv(df_cmc_t, paste0(x,"_cmc_aucell.csv"), row.names = TRUE)
}
