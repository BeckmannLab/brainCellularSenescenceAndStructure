ml R/4.4.1
R 

#libraries
library(dreamlet)
library(qvalue)
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)


#########
#up to 5
#########


var = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS"
), V2 = c("1percent", "5percent", "10percent", "20percent", "30percent", 
"1percent", "1percent", "5percent", "5percent", "10percent", 
"10percent", "20percent", "20percent", "30percent", "30percent", 
"1percent", "5percent", "10percent", "20percent", "30percent", 
"1percent", "1percent", "1percent", "5percent", "5percent", "5percent", 
"10percent", "10percent", "10percent", "20percent", "20percent", 
"20percent", "30percent", "30percent", "30percent"), V3 = c("pb_1percent", 
"pb_5percent", "pb_10percent", "pb_20percent", "pb_30percent", 
"pb_1percent", "pb_1percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_20percent", "pb_20percent", "pb_30percent", 
"pb_30percent", "pb_1percent", "pb_5percent", "pb_10percent", 
"pb_20percent", "pb_30percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent"
), V4 = c("Micro_TRUE", "Micro_TRUE", "Micro_TRUE", "Micro_TRUE", 
"Micro_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", 
"int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", 
"Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Oligo_TRUE", 
"OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", 
"Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", 
"Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE"), V5 = c("Micro_FALSE", 
"Micro_FALSE", "Micro_FALSE", "Micro_FALSE", "Micro_FALSE", "int_FALSE", 
"exc_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", 
"int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", "Vas_FALSE", 
"Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Oligo_FALSE", 
"OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", 
"Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", 
"Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE"), V6 = c("micro", 
"micro", "micro", "micro", "micro", "int", "exc", "int", "exc", 
"int", "exc", "int", "exc", "int", "exc", "vas", "vas", "vas", 
"vas", "vas", "Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", 
"Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC", 
"Astro"), V7 = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4
)), row.names = c(1L, 8L, 15L, 22L, 29L, 36L, 37L, 38L, 39L, 
40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 
53L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 62L, 63L, 64L, 65L
), class = "data.frame")

#pull the DE results into one object
all2 = c()
for (i in 1:nrow(var)) {
    file_path = paste0("correct_edger_upto5_", var[i,6], "_", var[i,2], "_fit.RDS")  
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        df$celltype = var[i,6]
        df$test = var[i,2]
        df$symbol = rownames(df)
        all2 = rbind(df, all2)
    }
}

sen_de = all2


##pull smri/snRNAseq results
var = structure(list(V1 = c("micro", "int", "exc", "vas", "Oligo", 
"OPC", "Astro"), V2 = c("30percent", "5percent", "5percent", 
"1percent", "10percent", "1percent", "1percent"), V3 = c("MG", 
"Int", "Exc", "NonNeu", "Oli", "OPC", "Ast")), class = "data.frame", row.names = c(NA, 
-7L))

data = #Supplementary Table 1 sn sMRI DE
myres <- c()
failures <- c()
for (i in unique(data$Feature)){
    name = i
    data <- data.table(feature=name, as.data.table(data))
        for (x in 1:nrow(var)){
            sen_de = all2[which(all2$celltype == var[x,1] & all2$test == var[x,2]),]
            sc_de = data[which(data$assay == var[x,3]),]
            p1_sc_de <- 1 - qvalue(sc_de$P.Value)$pi0
            p1_sen_de <- 1 - qvalue(sen_de$PValue)$pi0
             mer <- merge(sen_de, sc_de, by.x="symbol", by.y="ID")
             rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
             add <- data.table(feature=name, assay=var[x,1], pi1_sen_de=p1_sen_de,pi1_sc_de=p1_sc_de,spearmans=rho$estimate, p.value = rho$p.value)
            myres <- rbind(myres, add)
            print(add, row.names=FALSE, col.names="none")
            print("i", row.names=FALSE)    }
  } else {
    failures <- c(failures, i)
  }

saveRDS(myres, "cor_upto5_smri_aucell.RDS")

#############
###fiveandup
#############



var = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS", 
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", 
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", 
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", 
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", 
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS"
), V2 = c("1percent", "5percent", "10percent", "20percent", "30percent", 
"1percent", "1percent", "5percent", "5percent", "10percent", 
"10percent", "20percent", "20percent", "30percent", "30percent", 
"1percent", "5percent", "10percent", "20percent", "30percent", 
"1percent", "1percent", "1percent", "5percent", "5percent", "5percent", 
"10percent", "10percent", "10percent", "20percent", "20percent", 
"20percent", "30percent", "30percent", "30percent"), V3 = c("pb_1percent", 
"pb_5percent", "pb_10percent", "pb_20percent", "pb_30percent", 
"pb_1percent", "pb_1percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_20percent", "pb_20percent", "pb_30percent", 
"pb_30percent", "pb_1percent", "pb_5percent", "pb_10percent", 
"pb_20percent", "pb_30percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent"
), V4 = c("Micro_TRUE", "Micro_TRUE", "Micro_TRUE", "Micro_TRUE", 
"Micro_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", 
"int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", 
"Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE", "Oligo_TRUE", 
"OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", 
"Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", 
"Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE"), V5 = c("Micro_FALSE", 
"Micro_FALSE", "Micro_FALSE", "Micro_FALSE", "Micro_FALSE", "int_FALSE", 
"exc_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", 
"int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", "Vas_FALSE", 
"Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Oligo_FALSE", 
"OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", 
"Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", 
"Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE"), V6 = c("micro", 
"micro", "micro", "micro", "micro", "int", "exc", "int", "exc", 
"int", "exc", "int", "exc", "int", "exc", "vas", "vas", "vas", 
"vas", "vas", "Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", 
"Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC", 
"Astro"), V7 = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4
)), row.names = c(1L, 8L, 15L, 22L, 29L, 36L, 37L, 38L, 39L, 
40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 
53L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 62L, 63L, 64L, 65L
), class = "data.frame")

all2 = c()
for (i in 1:nrow(var2)){
    df = readRDS(paste0("correct_edger_5andup_",var2[i,6],"_",var2[i,2], "_fit.RDS"))
    df$celltype = var2[i,6]
    df$test = var2[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

sen_de = all2

##pull smri/snRNAseq results
var = structure(list(V1 = c("micro", "int", "exc", "vas", "Oligo", 
"OPC", "Astro"), V2 = c("30percent", "5percent", "5percent", 
"1percent", "10percent", "1percent", "1percent"), V3 = c("MG", 
"Int", "Exc", "NonNeu", "Oli", "OPC", "Ast")), class = "data.frame", row.names = c(NA, 
-7L))

data = #Supplementary Table 1 sn sMRI DE
myres <- c()
failures <- c()
for (i in unique(data$Feature)){
    name = i
    data <- data.table(feature=name, as.data.table(data))
        for (x in 1:nrow(var)){
            sen_de = all2[which(all2$celltype == var[x,1] & all2$test == var[x,2]),]
            sc_de = data[which(data$assay == var[x,3]),]
            p1_sc_de <- 1 - qvalue(sc_de$P.Value)$pi0
            p1_sen_de <- 1 - qvalue(sen_de$PValue)$pi0
             mer <- merge(sen_de, sc_de, by.x="symbol", by.y="ID")
             rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
             add <- data.table(feature=name, assay=var[x,1], pi1_sen_de=p1_sen_de,pi1_sc_de=p1_sc_de,spearmans=rho$estimate, p.value = rho$p.value)
            myres <- rbind(myres, add)
            print(add, row.names=FALSE, col.names="none")
            print("i", row.names=FALSE)    }
  } else {
    failures <- c(failures, i)
  }


saveRDS(myres, "cor_5andup_smri_aucell.RDS")

