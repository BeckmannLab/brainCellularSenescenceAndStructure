ml R/4.4.1
R
library(dreamlet)
library(qvalue)
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
library(limma)

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

all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("/correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

all2$sig = all2$FDR <= 0.05
sen_de = all2


##pull smri/snRNAseq results
var = structure(list(v1 = c("OPC", "Oligo", "Micro", "Int", "Exc", 
"Astro"), v2 = c("30percent", "30percent", "10percent", "30percent", 
"5percent", "30percent")), class = "data.frame", row.names = c(NA, 
-6L))

data = #Supplementary Table 1 sn sMRI DE
myres <- c()
failures <- c()
for (i in unique(data$Feature)){
    name = i
    data <- data.table(feature=name, as.data.table(data))
		for (x in 1:nrow(var)){
	  		sen_de = all2[which(all2$celltype == var[x,1] & all2$test == var[x,2]),]
	  		sc_de = data[which(data$new_celltype == var[x,1]),]
	  		p1_sc_de <- 1-propTrueNull(sc_de$P.Value)
	  		p1_sen_de <- 1-propTrueNull(sen_de$PValue)
					mer <- merge(sen_de, sc_de, by.x="symbol", by.y="ID")
					rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
					add <- data.table(feature=name, percent=var[x,2],celltype = var[x,1], pi1_sen_de=p1_sen_de,pi1_sc_de=p1_sc_de,spearmans=rho$estimate, p.value = rho$p.value)
	  		myres <- rbind(myres, add)
	  		print(add, row.names=FALSE, col.names="none")
	  		print("i", row.names=FALSE)    }
  } else {
    failures <- c(failures, i)
  }


saveRDS(myres, "corr_psychencodeAucellscore_smri_features.RDS")
