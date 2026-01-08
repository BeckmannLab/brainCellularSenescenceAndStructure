library(qs)
library(SingleCellExperiment)
library(Seurat)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)

sces = readRDS("sce_all_rosmap.RDS") #see setup

 unique(colData(sces)$broad.cell.type)
 "Ex"  "In"  "Opc" "Ast" "Oli" "End" "Mic" "Per"


colData(sces)$new_celltype <- ifelse(colData(sces)$broad.cell.type %in% c("Ex"), "Exc",
                               ifelse(colData(sces)$broad.cell.type %in% c("In"), "Int",
                               ifelse(colData(sces)$broad.cell.type %in% c("Mic"), "Micro",
                               ifelse(colData(sces)$broad.cell.type %in% c("Opc"), "OPC",
                               ifelse(colData(sces)$broad.cell.type %in% c("Oli"), "Oligo",
                               ifelse(colData(sces)$broad.cell.type %in% c("Ast"), "Astro",
                                      colData(sces)$broad.cell.type))))))

#subset by cell type 
celltypes = c("Exc", "Int", "Astro", "Oligo", "OPC", "Micro")
resfinal=foreach(x = celltypes) %dopar% {
	print(x)
	sub_cmc = sces[, which(colData(sces)$new_celltype == x)]
	df_cmc =counts(sub_cmc)
	df_cmc = as.matrix(df_cmc)
	df_cmc_t = t(df_cmc)
	write.csv(df_cmc_t, paste0("/sc/arion/projects/psychgen/lbp/data/neuroimaging/anina/brainSCOPE/expression/new_aucell/all_counts/by_celltype/",x,"_rosmap_aucell.csv"), row.names = TRUE)
}
