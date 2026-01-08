
##########################
##make all the bsub files 
##########################


for (x in 1:30) {
  i <- x
  
  # Create the content for each iteration, including the assignment of i
  r_content <- sprintf('
#ml R/4.4.1
#!/usr/bin/env Rscript

library(foreach)
library(dreamlet)
library(dplyr)
library(tidyr)
library(edgeR)
i = %d
folder <- "/final_pseudobulk/"
var <- structure(list(v1 = c("pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS", 
"pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS", 
"pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS", 
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS", 
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS", 
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS", 
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS", 
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS", 
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS", 
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS", 
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS", 
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS", 
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS", 
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS", 
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS"
), v2 = c("1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "30percent"), v3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_30percent", "pb_30percent"), v4 = c("astro_TRUE", 
"exc_TRUE", "int_TRUE", "micro_TRUE", "oli_TRUE", "opc_TRUE", 
"astro_TRUE", "exc_TRUE", "int_TRUE", "micro_TRUE", "oli_TRUE", 
"opc_TRUE", "astro_TRUE", "exc_TRUE", "int_TRUE", "micro_TRUE", 
"oli_TRUE", "opc_TRUE", "astro_TRUE", "exc_TRUE", "int_TRUE", 
"micro_TRUE", "oli_TRUE", "opc_TRUE", "astro_TRUE", "exc_TRUE", 
"int_TRUE", "micro_TRUE", "oli_TRUE", "opc_TRUE"), v5 = c("astro_FALSE", 
"exc_FALSE", "int_FALSE", "micro_FALSE", "oli_FALSE", "opc_FALSE", 
"astro_FALSE", "exc_FALSE", "int_FALSE", "micro_FALSE", "oli_FALSE", 
"opc_FALSE", "astro_FALSE", "exc_FALSE", "int_FALSE", "micro_FALSE", 
"oli_FALSE", "opc_FALSE", "astro_FALSE", "exc_FALSE", "int_FALSE", 
"micro_FALSE", "oli_FALSE", "opc_FALSE", "astro_FALSE", "exc_FALSE", 
"int_FALSE", "micro_FALSE", "oli_FALSE", "opc_FALSE"), v6 = c("Astro", 
"Exc", "Int", "Micro", "Oligo", "OPC", "Astro", "Exc", "Int", 
"Micro", "Oligo", "OPC", "Astro", "Exc", "Int", "Micro", "Oligo", 
"OPC", "Astro", "Exc", "Int", "Micro", "Oligo", "OPC", "Astro", 
"Exc", "Int", "Micro", "Oligo", "OPC"), v7 = c(0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4)), row.names = c(NA, -30L), class = "data.frame")
source("dreamletCompareClusters_edgeR.R")
source("processOneAssay_edgeR.R")
    print(paste(i, "start"))
    pb <- readRDS(paste0(folder,var[i,1]))
    ncells <- as.data.frame(int_colData(pb)$n_cells)
    ncells$celltype <- rownames(ncells) 

    ncells2 <- ncells %%>%%
      pivot_longer(!celltype,names_to = "id", values_to = "count")

    ncells2 <- as.data.frame(ncells2)

    ncells2$id <- gsub("\\\\.", "_", ncells2$id)

    colnames(ncells2) <- c("sen_auc", "Individual_ID", "ncell")

    metadata(pb)$aggr_means$Indv_ID <- metadata(pb)$aggr_means$Sample_ID 

    metadata(pb)$aggr_means$ncells <- as.numeric(ncells2$ncell)

    ct.pairs <- c(var[i,4], var[i,5])

    new_colnames <- make.names(colnames(assay(pb)))
    colnames(pb) <- new_colnames
    colnames(assay(pb)) <- new_colnames
    colData(pb)$Individual_ID = make.names(colData(pb)$Individual_ID)
    rownames(colData(pb)) = make.names(rownames(colData(pb)))
    metadata(pb)$aggr_means$Sample_ID = make.names(metadata(pb)$aggr_means$Sample_ID)
    metadata(pb)$aggr_means$Indv_ID = make.names(metadata(pb)$aggr_means$Indv_ID)


    try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    saveRDS(fit_top,paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))

} else {
    print(paste(i, "done"))
}', i)

  # Create the file name
  file_name <- sprintf("script_var%d.R", i)
  
  # Write the content to the file
  writeLines(r_content, file_name)
  
  print(paste("Created file:", file_name))
}

###########

for (x in 1:30) {
  r_content <- sprintf('#!/bin/bash
#BSUB -J R_job_var%d
#BSUB -o var%d_output.o
#BSUB -e var%d_output.e
#BSUB -q premium
#BSUB -P acc_psychgen
#BSUB -n 1
#BSUB -W 140:00 
#BSUB -R rusage[mem=8000]

ml R/4.4.1
Rscript script_var%d.R', x, x, x, x)

  # Create the file name
  file_name <- sprintf("script_var%d.sh", x)
  
  # Write the content to the file
  writeLines(r_content, file_name)
  
  print(paste("Created file:", file_name))
}

################
#submit jobs
cd /job_submissions/
for x in {1..30}
do
    bsub < "script_var${x}.sh"
done


####
#pull results

var = structure(list(v1 = c("pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS",
"pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS",
"pseudobulk_with_agg1percent.RDS", "pseudobulk_with_agg1percent.RDS",
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS",
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS",
"pseudobulk_with_agg5percent.RDS", "pseudobulk_with_agg5percent.RDS",
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS",
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS",
"pseudobulk_with_agg10percent.RDS", "pseudobulk_with_agg10percent.RDS",
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS",
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS",
"pseudobulk_with_agg20percent.RDS", "pseudobulk_with_agg20percent.RDS",
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS",
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS",
"pseudobulk_with_agg30percent.RDS", "pseudobulk_with_agg30percent.RDS"
), v2 = c("1percent", "1percent", "1percent", "1percent", "1percent",
"1percent", "5percent", "5percent", "5percent", "5percent", "5percent",
"5percent", "10percent", "10percent", "10percent", "10percent",
"10percent", "10percent", "20percent", "20percent", "20percent",
"20percent", "20percent", "20percent", "30percent", "30percent",
"30percent", "30percent", "30percent", "30percent"), v3 = c("pb_1percent",
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent",
"pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent",
"pb_5percent", "pb_10percent", "pb_10percent", "pb_10percent",
"pb_10percent", "pb_10percent", "pb_10percent", "pb_20percent",
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_30percent", "pb_30percent", "pb_30percent"), v4 = c("astro_TRUE",
"exc_TRUE", "int_TRUE", "micro_TRUE", "oli_TRUE", "opc_TRUE",
"astro_TRUE", "exc_TRUE", "int_TRUE", "micro_TRUE", "oli_TRUE",
"opc_TRUE", "astro_TRUE", "exc_TRUE", "int_TRUE", "micro_TRUE",
"oli_TRUE", "opc_TRUE", "astro_TRUE", "exc_TRUE", "int_TRUE",
"micro_TRUE", "oli_TRUE", "opc_TRUE", "astro_TRUE", "exc_TRUE",
"int_TRUE", "micro_TRUE", "oli_TRUE", "opc_TRUE"), v5 = c("astro_FALSE",
"exc_FALSE", "int_FALSE", "micro_FALSE", "oli_FALSE", "opc_FALSE",
"astro_FALSE", "exc_FALSE", "int_FALSE", "micro_FALSE", "oli_FALSE",
"opc_FALSE", "astro_FALSE", "exc_FALSE", "int_FALSE", "micro_FALSE",
"oli_FALSE", "opc_FALSE", "astro_FALSE", "exc_FALSE", "int_FALSE",
"micro_FALSE", "oli_FALSE", "opc_FALSE", "astro_FALSE", "exc_FALSE",
"int_FALSE", "micro_FALSE", "oli_FALSE", "opc_FALSE"), v6 = c("Astro",
"Exc", "Int", "Micro", "Oligo", "OPC", "Astro", "Exc", "Int",
"Micro", "Oligo", "OPC", "Astro", "Exc", "Int", "Micro", "Oligo",
"OPC", "Astro", "Exc", "Int", "Micro", "Oligo", "OPC", "Astro",
"Exc", "Int", "Micro", "Oligo", "OPC"), v7 = c(0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4)), row.names = c(NA, -30L), class = "data.frame")
all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}
