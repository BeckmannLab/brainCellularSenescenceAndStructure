
ml R/4.4.1
R
library(dreamlet)
library(qvalue)
library(data.table)

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


all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

liv_aucell = all2

liv_aucell$celltype = gsub("noneu","nonneu", liv_aucell$celltype)

var = structure(list(V1 = c("oli", "oli", "oli", "oli", "oli", "oli", 
"oli", "exc", "exc", "exc", "exc", "exc", "exc", "exc", "int", 
"int", "int", "int", "int", "int", "int", "nonneu", "nonneu", 
"nonneu", "nonneu", "nonneu", "nonneu", "nonneu", "mg", "mg", 
"mg", "mg", "mg", "mg", "mg", "ast", "ast", "ast", "ast", "ast", 
"ast", "ast", "opc", "opc", "opc", "opc", "opc", "opc", "opc"
), V2 = c("10percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "5percent", "5percent", "5percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "1percent", "1percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "1percent", "1percent", "1percent", "1percent"
), V3 = c("Oli", "Oli", "Oli", "Oli", "Oli", "Oli", "Oli", "Exc", 
"Exc", "Exc", "Exc", "Exc", "Exc", "Exc", "Int", "Int", "Int", 
"Int", "Int", "Int", "Int", "NonNeu", "NonNeu", "NonNeu", "NonNeu", 
"NonNeu", "NonNeu", "NonNeu", "MG", "MG", "MG", "MG", "MG", "MG", 
"MG", "Ast", "Ast", "Ast", "Ast", "Ast", "Ast", "Ast", "OPC", 
"OPC", "OPC", "OPC", "OPC", "OPC", "OPC"), V4 = c("oli", "exc", 
"int", "nonneu", "mg", "ast", "opc", "oli", "exc", "int", "nonneu", 
"mg", "ast", "opc", "oli", "exc", "int", "nonneu", "mg", "ast", 
"opc", "oli", "exc", "int", "nonneu", "mg", "ast", "opc", "oli", 
"exc", "int", "nonneu", "mg", "ast", "opc", "oli", "exc", "int", 
"nonneu", "mg", "ast", "opc", "oli", "exc", "int", "nonneu", 
"mg", "ast", "opc"), V5 = c("10percent", "5percent", "5percent", 
"1percent", "30percent", "1percent", "1percent", "10percent", 
"5percent", "5percent", "1percent", "30percent", "1percent", 
"1percent", "10percent", "5percent", "5percent", "1percent", 
"30percent", "1percent", "1percent", "10percent", "5percent", 
"5percent", "1percent", "30percent", "1percent", "1percent", 
"10percent", "5percent", "5percent", "1percent", "30percent", 
"1percent", "1percent", "10percent", "5percent", "5percent", 
"1percent", "30percent", "1percent", "1percent", "10percent", 
"5percent", "5percent", "1percent", "30percent", "1percent", 
"1percent"), V6 = c("Oli", "Exc", "Int", "NonNeu", "MG", "Ast", 
"OPC", "Oli", "Exc", "Int", "NonNeu", "MG", "Ast", "OPC", "Oli", 
"Exc", "Int", "NonNeu", "MG", "Ast", "OPC", "Oli", "Exc", "Int", 
"NonNeu", "MG", "Ast", "OPC", "Oli", "Exc", "Int", "NonNeu", 
"MG", "Ast", "OPC", "Oli", "Exc", "Int", "NonNeu", "MG", "Ast", 
"OPC", "Oli", "Exc", "Int", "NonNeu", "MG", "Ast", "OPC")), class = "data.frame", row.names = c(NA, 
-49L))


myres <- c()
for (i in 1:nrow(var)){
  cell_1 = liv_aucell[liv_aucell$celltype == var[i,1] &liv_aucell$test == var[i,2],]
  cell_2 = liv_aucell[liv_aucell$celltype == var[i,4] & liv_aucell$test == var[i,5],]
    mer <- merge(cell_1, cell_2, by="symbol")
  rho <- cor.test(mer$logFC.x, mer$logFC.y, method="spearman")
  add <- data.table(celltype1 = var[i,1], celltype2 = var[i,4], spearmans=rho$estimate, p.value = rho$p.value)
  myres <- rbind(myres, add)
  print(i)
}

liv_cor = myres
liv_cor$ds = "liv"

liv_cor = liv_cor[liv_cor$celltype1 != "nonneu",]
liv_cor = liv_cor[liv_cor$celltype2 != "nonneu",]

liv_cor$celltype1 = gsub("mg", "MG", liv_cor$celltype1)
liv_cor$celltype1 = gsub("opc", "OPC", liv_cor$celltype1)
liv_cor$celltype1 = gsub("exc", "Exc", liv_cor$celltype1)
liv_cor$celltype1 = gsub("int", "Inh", liv_cor$celltype1)
liv_cor$celltype1 = gsub("oli", "Oli", liv_cor$celltype1)
liv_cor$celltype1 = gsub("ast", "Ast", liv_cor$celltype1)

liv_cor$celltype2 = gsub("mg", "MG", liv_cor$celltype2)
liv_cor$celltype2 = gsub("opc", "OPC", liv_cor$celltype2)
liv_cor$celltype2 = gsub("exc", "Exc", liv_cor$celltype2)
liv_cor$celltype2 = gsub("int", "Inh", liv_cor$celltype2)
liv_cor$celltype2 = gsub("oli", "Oli", liv_cor$celltype2)
liv_cor$celltype2 = gsub("ast", "Ast", liv_cor$celltype2)


g = ggplot(liv_cor, aes(x = celltype1, y = celltype2, fill = spearmans)) +
  geom_tile() +
  geom_text(aes(label = paste0(
    ifelse(p.value < 2.2e-308, "< 2.2e-308", formatC(p.value, format = "e", digits = 2))
  )), size = 6) +  # Display "< 2.2e-308" for very small p-values, otherwise in scientific notation
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", 
    midpoint = 0, limit = c(0, 1), space = "Lab", 
    name = expression("Spearman's rho (" * rho * ")")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 18, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.key.size = unit(1, "cm")
  ) + 
  labs(y = "Cell Type 1", x = "Cell Type 2")
ggsave("cor_across_cells_liv.pdf", g, width = 11,height = 9)
