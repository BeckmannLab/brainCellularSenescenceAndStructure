#combine the celtypes
library(data.table)
library(doParallel)
registerDoParallel(cores=40)

#see aucell_by_celltype
micro = fread("micro_aucell_across.csv")
astro = fread("astro_aucell_across.csv")
exc = fread("exc_aucell_across.csv")
int = fread("int_aucell_across.csv")
opc = fread("OPC_aucell_across.csv")
oli = fread("oli_aucell_across.csv")

celltypes = c("micro", "astro", "exc", "int", "opc", "oli")

all_var = structure(list(V1 = c("1percent", "5percent", "10percent", "20percent", 
"30percent"), V2 = c(0.99, 0.95, 0.9, 0.8, 0.7)), row.names = c(NA, 
-5L), class = "data.frame")

foreach(x = celltypes)%dopar%{
	print(x)
	final_df <- data.frame()
	for (i in 1:nrow(all_var)){
        AUC_mat <- get(x)
        threshold <- quantile(AUC_mat$MyGeneSet, all_var[i,2])
        AUC_mat$sen <- AUC_mat$MyGeneSet >= threshold
        AUC_mat$percent = all_var[i,1]
        AUC_mat$celltype = x
        print(i)
        final_df = rbind(AUC_mat, final_df)
    }
	saveRDS(final_df,paste0("thresholds_across_samples_", x, ".RDS"))
}

##check 
exc_aucell = readRDS("thresholds_across_samples_exc.RDS")
int_aucell = readRDS("thresholds_across_samples_int.RDS")
micro_aucell = readRDS("thresholds_across_samples_micro.RDS")
oli_aucell = readRDS("thresholds_across_samples_oli.RDS")
opc_aucell = readRDS("thresholds_across_samples_opc.RDS")
astro_aucell = readRDS("thresholds_across_samples_astro.RDS")

all_sene = rbind(exc_aucell, int_aucell, micro_aucell, oli_aucell, opc_aucell, astro_aucell)
saveRDS(all_sene, "across_all_celltypes_and_thresholds.RDS")
