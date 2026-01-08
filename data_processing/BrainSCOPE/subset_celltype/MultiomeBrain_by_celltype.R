library(qs)
library(SingleCellExperiment)
library(Seurat)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)

sces = qread("sce_all_MultiomeBrain.qs") #see setup

#subset by cell type 
celltypes = c("Exc", "Int", "Astro", "Oligo", "OPC", "Micro")
resfinal=foreach(x = celltypes) %dopar% {
	print(x)
	sub_cmc = sces[, which(colData(sces)$new_celltype == x)]
	df_cmc =counts(sub_cmc)
	df_cmc = as.matrix(df_cmc)
	df_cmc_t = t(df_cmc)
	write.csv(df_cmc_t, paste0(x,"_MultiomeBrain_aucell.csv"), row.names = TRUE)
}
