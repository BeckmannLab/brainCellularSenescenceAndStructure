###og
library(qvalue)
library(data.table)
library(seriation)
library(corrplot)

var = structure(list(V1 = c("lbp_06.21.24_pseudobulk_1percent.RDS",
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS",
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS",
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS",
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS",
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS",
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS",
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS",
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS",
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS",
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS",
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS",
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS",
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS",
"lbp_06.21.24_pseudobulk_30percent.RDS"), V2 = c("1percent",
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent",
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent",
"5percent", "20percent", "20percent", "20percent", "20percent",
"20percent", "20percent", "20percent", "30percent", "30percent",
"30percent", "30percent", "30percent", "30percent", "30percent"
), V3 = c("pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent",
"pb_1percent", "pb_1percent", "pb_1percent", "pb_5percent", "pb_5percent",
"pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent",
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_20percent", "pb_20percent", "pb_30percent",
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_30percent", "pb_30percent"), V4 = c("Ast_TRUE", "Exc_TRUE",
"Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE",
"Ast_FALSE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE",
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE",
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE", "Ast_FALSE",
"Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE",
"OPC_TRUE"), V5 = c("Ast_FALSE", "Exc_FALSE", "Int_FALSE", "MG_FALSE",
"NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE",
"Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE",
"Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE",
"Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE",
"MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE"), V6 = c("ast",
"exc", "int", "mg", "noneu", "oli", "opc", "ast", "exc", "int",
"mg", "noneu", "oli", "opc", "ast", "exc", "int", "mg", "noneu",
"oli", "opc", "ast", "exc", "int", "mg", "noneu", "oli", "opc"
), V7 = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4)), class = "data.frame", row.names = c(NA,
-28L))

var2 = structure(list(V1 = c("lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS"
), V2 = c("10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent"), V3 = c("pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent"), V4 = c("Ast_FALSE", "Exc_TRUE", 
"Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE"), 
    V5 = c("Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", 
    "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE"), V6 = c("ast", 
    "exc", "int", "mg", "noneu", "oli", "opc"), V7 = c(0.4, 0.4, 
    0.4, 0.4, 0.4, 0.4, 0.4)), class = "data.frame", row.names = c(NA, 
-7L))
var = rbind(var, var2)
var_for_og = var[c(1,9,10,25,34,28),]

all2 = c()
for (i in 1:nrow(var_for_og)){
    df = readRDS(paste0("correct_edger_",var_for_og[i,6],"_",var_for_og[i,2], "_fit.RDS"))
    df$celltype = var_for_og[i,6]
    df$test = var_for_og[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

liv_aucell = all2


#nd up to five

var = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS",
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS",
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS",
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS",
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS",
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS",
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS",
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
), V2 = c("1percent", "1percent", "1percent", "1percent", "1percent",
"1percent", "1percent", "5percent", "5percent", "5percent", "5percent",
"5percent", "5percent", "5percent", "10percent", "10percent",
"10percent", "10percent", "10percent", "10percent", "10percent",
"20percent", "20percent", "20percent", "20percent", "20percent",
"20percent", "20percent", "30percent", "30percent", "30percent",
"30percent", "30percent", "30percent", "30percent", "1percent",
"1percent", "5percent", "5percent", "10percent", "10percent",
"20percent", "20percent", "30percent", "30percent", "1percent",
"5percent", "10percent", "20percent", "30percent", "1percent",
"1percent", "1percent", "5percent", "5percent", "5percent", "10percent", 
"10percent", "10percent", "20percent", "20percent", "20percent", 
"30percent", "30percent", "30percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_1percent", "pb_1percent", "pb_5percent", "pb_5percent", "pb_10percent",
"pb_10percent", "pb_20percent", "pb_20percent", "pb_30percent",
"pb_30percent", "pb_1percent", "pb_5percent", "pb_10percent",
"pb_20percent", "pb_30percent", "pb_1percent", "pb_1percent",
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent",
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent"), V4 = c("Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE",
"L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE", "Micro_TRUE",
"L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE",
"LAMP5_NOS1_TRUE", "PN_dev_TRUE", "Micro_TRUE", "L2-3_CUX2_TRUE",
"L4_RORB_TRUE", "L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE",
"PN_dev_TRUE", "Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE",
"L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE",
"Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE",
"L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE", "int_TRUE",
"exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE",
"exc_TRUE", "int_TRUE", "exc_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE",
"Vas_TRUE", "Vas_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE",
"Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE",
"Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE",
"OPC_TRUE", "Astro_TRUE"), V5 = c("Micro_FALSE", "L2-3_CUX2_FALSE",
"L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE",
"PN_dev_FALSE", "Micro_FALSE", "L2-3_CUX2_FALSE", "L4_RORB_FALSE",
"L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE", "PN_dev_FALSE",
"Micro_FALSE", "L2-3_CUX2_FALSE", "L4_RORB_FALSE", "L5-6_THEMIS_FALSE",
"L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE", "PN_dev_FALSE", "Micro_FALSE",
"L2-3_CUX2_FALSE", "L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE",
"LAMP5_NOS1_FALSE", "PN_dev_FALSE", "Micro_FALSE", "L2-3_CUX2_FALSE",
"L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE",
"PN_dev_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE",
"int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE",
"exc_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE",
"Vas_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE",
"OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE",
"Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE",
"Astro_FALSE"), V6 = c("micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS",
"L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2", "L4_RORB",
"L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2",
"L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "micro",
"L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5",
"PN_dev", "micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4",
"LAMP5", "PN_dev", "int", "exc", "int", "exc", "int", "exc",
"int", "exc", "int", "exc", "vas", "vas", "vas", "vas", "vas",
"Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC",
"Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro"), V7 = c(0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)), class = "data.frame", row.names = c(NA,
-65L))

var_upto5 = var[c(65,41,36,22,63,55),]


all2 = c()
for (i in 1:nrow(var_upto5)) {
    file_path = paste0("correct_edger_upto5_", var_upto5[i,6], "_", var_upto5[i,2], "_fit.RDS")  
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        df$celltype = var_upto5[i,6]
        df$test = var_upto5[i,2]
        df$symbol = rownames(df)
        all2 = rbind(df, all2)
    }
}
nd_upto5_aucell = all2


#nd five and up

var = structure(list(V1 = c("pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS",
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS",
"pseudobulk_1percent.RDS", "pseudobulk_1percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS", "pseudobulk_5percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS", "pseudobulk_10percent.RDS",
"pseudobulk_10percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS",
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS",
"pseudobulk_20percent.RDS", "pseudobulk_20percent.RDS", "pseudobulk_30percent.RDS",
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS",
"pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS", "pseudobulk_30percent.RDS",
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
), V2 = c("1percent", "1percent", "1percent", "1percent", "1percent",
"1percent", "1percent", "5percent", "5percent", "5percent", "5percent",
"5percent", "5percent", "5percent", "10percent", "10percent",
"10percent", "10percent", "10percent", "10percent", "10percent",
"20percent", "20percent", "20percent", "20percent", "20percent",
"20percent", "20percent", "30percent", "30percent", "30percent",
"30percent", "30percent", "30percent", "30percent", "1percent",
"1percent", "5percent", "5percent", "10percent", "10percent",
"20percent", "20percent", "30percent", "30percent", "1percent",
"5percent", "10percent", "20percent", "30percent", "1percent",
"1percent", "1percent", "5percent", "5percent", "5percent", "10percent", 
"10percent", "10percent", "20percent", "20percent", "20percent", 
"30percent", "30percent", "30percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent",
"pb_1percent", "pb_1percent", "pb_5percent", "pb_5percent", "pb_10percent",
"pb_10percent", "pb_20percent", "pb_20percent", "pb_30percent",
"pb_30percent", "pb_1percent", "pb_5percent", "pb_10percent",
"pb_20percent", "pb_30percent", "pb_1percent", "pb_1percent",
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_10percent",
"pb_10percent", "pb_10percent", "pb_20percent", "pb_20percent",
"pb_20percent", "pb_30percent", "pb_30percent", "pb_30percent"
),V4 = c("Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE",
"L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE", "Micro_TRUE",
"L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE",
"LAMP5_NOS1_TRUE", "PN_dev_TRUE", "Micro_TRUE", "L2-3_CUX2_TRUE",
"L4_RORB_TRUE", "L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE",
"PN_dev_TRUE", "Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE",
"L5-6_THEMIS_TRUE", "L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE",
"Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE",
"L5-6_TLE4_TRUE", "LAMP5_NOS1_TRUE", "PN_dev_TRUE", "int_TRUE",
"exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE", "exc_TRUE", "int_TRUE",
"exc_TRUE", "int_TRUE", "exc_TRUE", "Vas_TRUE", "Vas_TRUE", "Vas_TRUE",
"Vas_TRUE", "Vas_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE",
"Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE", "OPC_TRUE",
"Astro_TRUE", "Oligo_TRUE", "OPC_TRUE", "Astro_TRUE", "Oligo_TRUE",
"OPC_TRUE", "Astro_TRUE"), V5 = c("Micro_FALSE", "L2-3_CUX2_FALSE",
"L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE",
"PN_dev_FALSE", "Micro_FALSE", "L2-3_CUX2_FALSE", "L4_RORB_FALSE",
"L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE", "PN_dev_FALSE",
"Micro_FALSE", "L2-3_CUX2_FALSE", "L4_RORB_FALSE", "L5-6_THEMIS_FALSE",
"L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE", "PN_dev_FALSE", "Micro_FALSE",
"L2-3_CUX2_FALSE", "L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE",
"LAMP5_NOS1_FALSE", "PN_dev_FALSE", "Micro_FALSE", "L2-3_CUX2_FALSE",
"L4_RORB_FALSE", "L5-6_THEMIS_FALSE", "L5-6_TLE4_FALSE", "LAMP5_NOS1_FALSE",
"PN_dev_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE",
"int_FALSE", "exc_FALSE", "int_FALSE", "exc_FALSE", "int_FALSE",
"exc_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE", "Vas_FALSE",
"Vas_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE",
"OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE", "Astro_FALSE",
"Oligo_FALSE", "OPC_FALSE", "Astro_FALSE", "Oligo_FALSE", "OPC_FALSE",
"Astro_FALSE"), V6 = c("micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS",
"L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2", "L4_RORB",
"L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2",
"L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "micro",
"L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5",
"PN_dev", "micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4",
"LAMP5", "PN_dev", "int", "exc", "int", "exc", "int", "exc",
"int", "exc", "int", "exc", "vas", "vas", "vas", "vas", "vas",
"Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC",
"Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro"), V7 = c(0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)), class = "data.frame", row.names = c(NA,
-65L))
var_5and_up = var[c(62,41,36,29,54,64),]

all2 = c()
for (i in 1:nrow(var_5and_up)){
    df = readRDS(paste0("correct_edger_5andup_",var_5and_up[i,6],"_",var_5and_up[i,2], "_fit.RDS"))
    df$celltype = var_5and_up[i,6]
    df$test = var_5and_up[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

nd_5andup_aucell = all2


##correlations

#og to nd up to five
var_for_og2 = var_for_og[,c(2,6)]
var_upto52 = var_upto5[,c(2,6)]
colnames(var_upto52) = c("v2", "v6")
var_for_og2 = var_for_og2[order(var_for_og2$V6),]
var_upto52 = var_upto52[order(var_upto52$v6),]

all_var = cbind(var_upto52,var_for_og2)

myres <- c()

for (x in 1:nrow(all_var)){
    liv_aucell_subset = liv_aucell[which(liv_aucell$celltype == all_var[x,4] & liv_aucell$test == all_var[x,3]),]
    nd_upto5_aucell_subset = nd_upto5_aucell[which(nd_upto5_aucell$celltype == all_var[x,2] & nd_upto5_aucell$test == all_var[x,1]),]
    mer <- merge(liv_aucell_subset, nd_upto5_aucell_subset, by="symbol")
    rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
    p_val = rho$p.value
	fdr = p.adjust(p_val, method = "fdr")
    add <- data.table(celltype =all_var[x,2], spearmans=rho$estimate, p = p_val, fdr = fdr )
        myres <- rbind(myres, add)
        print(add, row.names=FALSE, col.names="none") 
    } 

og_nd_upto5 = myres
og_nd_upto5$test = "og_nd_upto5"

#og to nd five and up
var_for_og2 = var_for_og[,c(2,6)]
var_5and_up2 = var_5and_up[,c(2,6)]
colnames(var_5and_up2) = c("v2", "v6")
var_for_og2 = var_for_og2[order(var_for_og2$V6),]
var_5and_up2 = var_5and_up2[order(var_5and_up2$v6),]

all_var = cbind(var_5and_up2,var_for_og2)

myres <- c()

for (x in 1:nrow(all_var)){
    liv_aucell_subset = liv_aucell[which(liv_aucell$celltype == all_var[x,4] & liv_aucell$test == all_var[x,3]),]
    nd_5andup_aucell_subset = nd_5andup_aucell[which(nd_5andup_aucell$celltype == all_var[x,2] & nd_5andup_aucell$test == all_var[x,1]),]
    mer <- merge(liv_aucell_subset, nd_5andup_aucell_subset, by="symbol")
    rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
    p_val = rho$p.value
	fdr = p.adjust(p_val, method = "fdr")
    add <- data.table(celltype =all_var[x,2], spearmans=rho$estimate, p = p_val, fdr = fdr )
        myres <- rbind(myres, add)
        print(add, row.names=FALSE, col.names="none") 
    } 

og_nd_5andup = myres
og_nd_5andup$test = "og_nd_5andup"

###plot

all_cors = rbind(og_nd_upto5,og_nd_5andup)
all_cors$new_celltype = all_cors$celltype
all_cors$new_celltype = gsub("Astro", "Astrocytes",all_cors$new_celltype )
all_cors$new_celltype = gsub("ast", "Astrocytes",all_cors$new_celltype )
all_cors$new_celltype = gsub("exc", "Excitatory Neurons",all_cors$new_celltype )
all_cors$new_celltype = gsub("Exc", "Excitatory Neurons",all_cors$new_celltype )
all_cors$new_celltype = gsub("int", "Inhibitory Neurons",all_cors$new_celltype )
all_cors$new_celltype = gsub("Int", "Inhibitory Neurons",all_cors$new_celltype )
all_cors$new_celltype = gsub("micro", "Microglia",all_cors$new_celltype )
all_cors$new_celltype = gsub("Micro", "Microglia",all_cors$new_celltype )
all_cors$new_celltype = gsub("Oligo", "Oligodendrocytes",all_cors$new_celltype )
all_cors$new_celltype = gsub("oli", "Oligodendrocytes",all_cors$new_celltype )
all_cors$new_celltype = gsub("OPC", "Oligodendrocyte\nProgenitor Cells",all_cors$new_celltype )
all_cors$new_celltype = gsub("opc", "Oligodendrocyte\nProgenitor Cells",all_cors$new_celltype )

all_cors$new_celltype = gsub("Excitatory Neuronsitatory Neurons", "Excitatory Neurons",all_cors$new_celltype )
all_cors$new_celltype = gsub("Microgliaglia", "Microglia",all_cors$new_celltype )
all_cors$new_celltype = gsub("mg", "Microglia",all_cors$new_celltype )

all_cors$test = gsub("og_nd_upto5","ND up to five",all_cors$test)
all_cors$test = gsub("og_nd_5andup","ND five and up",all_cors$test)


library(ggplot2)
g = ggplot(all_cors, aes(x = reorder(new_celltype, -spearmans), 
                         y = spearmans, 
                         color = test)) +
  geom_point(size = 6) +
  theme_bw() +
  labs(x = "Cell Type", 
       y = "Correlation (Spearman's \u03C1)",
       color = "Dataset") +
  # Use alternative colorblind-friendly palette
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  theme(
    # Position the legend inside the top-right corner
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    # Add border (box) around the legend
    legend.box.background = element_rect(colour = "black"),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 17)
  )

ggsave("nd_sen_correlation.pdf",g, width = 9, height = 7)
