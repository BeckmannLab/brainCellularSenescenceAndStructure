library(qs)
library(SingleCellExperiment)
library(Seurat)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)

sces = readRDS("sce_all_SZBDMulti-Seq.RDS") #see setup

#subset by cell type 
celltypes = c("Exc", "Int", "Astro", "Oligo", "OPC", "Micro")
resfinal=foreach(i = names(sces)) %dopar% {
    foreach(x = celltypes)%dopar% {
    print(x)
    print(i)
    name = sub("-.*", "", i)
    sces3 <- sces[[i]]
    colData(sces3)$new_celltype <- ifelse(colData(sces3)$cellType %in% 
                                    c("L2/3 IT", "L4 IT", "L5 IT", "L5 ET", 
                                      "L5/6 NP", "L6 IT", "L6 IT Car3", 
                                      "L6 CT", "L6b"), "Exc", 
                                    ifelse(colData(sces3)$cellType %in% 
                                           c("Pvalb", "Sst", "Sst Chodl", "Lamp5", 
                                             "Lamp5 Lhx6", "Chandelier", "Vip", 
                                             "Sncg", "Pax6"), "Int", 
                                           colData(sces3)$cellType))
    sub_cmc = sces3[, which(colData(sces3)$new_celltype == x)]
    df_cmc =counts(sub_cmc)
    df_cmc = as.matrix(df_cmc)
    df_cmc_t = t(df_cmc)
    write.csv(df_cmc_t, paste0(x,"_", name, "_SZBDMulti-Seq_aucell.csv"), row.names = TRUE)
}}
