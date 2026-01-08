library(data.table)
library(ggpubr)
library(writexl)
library(dplyr)

senCID = readRDS("all_sc_senCID_output.RDS") #senCID output

mg = #senCID results 
exc = #senCID results 
int = #senCID results 
opc = #senCID results 
ast = #senCID results 
oli = #senCID results 

aucell_results = rbind(mg, exc, int, nonneu, opc, ast, oli)
aucell_results$V1 = rownames(aucell_results)

final_all = merge(aucell_results, senCID, by = "V1")

dirs <- list.dirs(path = "path", full.names = TRUE, recursive = FALSE)

dirs = dirs[-1]

df_final=c()
for (x in dirs){
	print(x)
	 # files <- list.files(dirs)
	 resSID = fread(paste0(x, "/recSID.csv"))
	 df_final = rbind(resSID,df_final)
}

aucell = aucell_results[,c("celltype", "percent")]
aucell$V1 = rownames(aucell)

all = merge(df_final, aucell, by = "V1")
all = as.data.frame(all)
all2 = merge(all[,c(1,2)], final_all, by.x = c("V1","RecSID"), by.y = c("V1", "SID"))

correlations <- all2 %>%
  group_by(celltype) %>%
  summarise(
    cor_test = list(cor.test(Decision, sen_list, method = "spearman")),
    rho = cor_test[[1]]$estimate[[1]],  # Extract rho properly
    p_value = cor_test[[1]]$p.value  # Extract p-value properly
  ) %>%
  ungroup() %>%
  mutate(
    fdr = p.adjust(p_value, method = "fdr"),
    fdr_label = formatC(fdr, format = "e", digits = 2),  # Use scientific notation
    label = sprintf("rho = %.2f, FDR = %s", rho, fdr_label)  # Concatenate rho and FDR
  )


correlations$new_celltype = correlations$celltype 
correlations$new_celltype = gsub("Exc", "Excitatory Neurons",correlations$new_celltype )
correlations$new_celltype = gsub("Int", "Inhibitory Neurons",correlations$new_celltype )
correlations$new_celltype = gsub("MG", "Microglia",correlations$new_celltype )
correlations$new_celltype = gsub("Ast", "Astrocytes",correlations$new_celltype )
correlations$new_celltype = gsub("Oli", "Oligodendrocytes",correlations$new_celltype )
correlations$new_celltype = gsub("OPC", "Oligodendrocyte progenitor cells",correlations$new_celltype )

all3$new_celltype = all3$celltype 
all3$new_celltype = gsub("Exc", "Excitatory Neurons",all3$new_celltype )
all3$new_celltype = gsub("Int", "Inhibitory Neurons",all3$new_celltype )
all3$new_celltype = gsub("MG", "Microglia",all3$new_celltype )
all3$new_celltype = gsub("Ast", "Astrocytes",all3$new_celltype )
all3$new_celltype = gsub("Oli", "Oligodendrocytes",all3$new_celltype )
all3$new_celltype = gsub("OPC", "Oligodendrocyte progenitor cells",all3$new_celltype )

all4 = all3[,c("V1", "sen_list", "new_celltype", "Decision")]

write_xlsx(all4, "sencid_aucell_comaprison.xlsx")
