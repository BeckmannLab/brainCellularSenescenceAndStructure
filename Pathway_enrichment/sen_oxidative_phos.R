
#oxidative phosphorylation
wget https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_OXIDATIVE_PHOSPHORYLATION&fileType=gmt


kegg = read.delim("hallmark_oxidative_phosphorylation_msigdb_terms.gmt")
kegg = as.data.frame(kegg[-1,])
colnames(kegg) = "symbol"

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
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS"
), V2 = c("1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent"
), V4 = c("Ast_TRUE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", 
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE", "Ast_FALSE", 
"Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", 
"OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", 
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE"), V5 = c("Ast_FALSE", 
"Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", 
"OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", 
"NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", 
"Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", 
"Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", 
"Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE", 
"MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE"), V6 = c("ast", 
"exc", "int", "mg", "noneu", "oli", "opc", "ast", "exc", "int", 
"mg", "noneu", "oli", "opc", "ast", "exc", "int", "mg", "noneu", 
"oli", "opc", "ast", "exc", "int", "mg", "noneu", "oli", "opc", 
"ast", "exc", "int", "mg", "noneu", "oli", "opc"), V7 = c(0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)), row.names = c(NA, -35L
), class = "data.frame")

all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("correct_edger_",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}


all2$sig_logfc_less0 <- (all2$FDR <= 0.05) & (all2$logFC < 0)
all2$sig_logfc_greater0 <- (all2$FDR <= 0.05) & (all2$logFC > 0)
kegg$status = "kegg_g"
sc_with_kegg=merge(all2,kegg, by = "symbol", all.x = T)
sc_with_kegg$status[is.na(sc_with_kegg$status)]="not_kegg_g"


var = structure(list(V1 = c("1percent", "10percent", "5percent", "5percent", 
"30percent", "30percent"), V2 = c("ast", "oli", "exc", "int", 
"mg", "opc")), class = "data.frame", row.names = c(NA, -6L))

sc_with_kegg$status <- factor(sc_with_kegg$status, c("not_kegg_g", "kegg_g"))

#test
results_greater_logfcgreater0 = list() 
for (i in 1:nrow(var)) {
    sc_subset = sc_with_kegg[which(sc_with_kegg$celltype == var[i,2] & sc_with_kegg$test == var[i,1]), ]
    test = tryCatch(fisher.test(table(sc_subset$status, sc_subset$sig_logfc_greater0), alternative = "greater"), error = function(e) list())  
    # Check if the test is not an empty list (i.e., the test was successful)
    if (length(test) > 0) {
        results_greater_logfcgreater0[[i]] = list(
            var1 = var[i,1], 
            var2 = var[i,2], 
            table = table(sc_subset$status, sc_subset$sig_logfc_greater0), 
            test = test, 
            total_expressed = length(sc_subset$celltype), 
            degs = length(which(sc_subset$FDR <= 0.05)), 
            kegg_TRUE = length(which(sc_subset$status == "kegg_g" & sc_subset$sig_logfc_greater0 == TRUE)), 
            total_expressed_kegg = length(which(sc_subset$status == "kegg_g"))
        )
    }
}

final_df_greater_logfcgreater0 = data.frame(var1 = character(),
                      var2 = character(),
                      p_value = numeric(),
                      odds_ratio= numeric(),
                      total_expressed = numeric(),
                      degs = numeric(),
                      kegg_true = numeric(),
                      total_expressed_kegg = numeric(),
                      fdr = numeric())

# Iterate through the results list and extract the necessary information
for (i in 1:length(results_greater_logfcgreater0)) {
    result = results_greater_logfcgreater0[[i]]
    
    # Extract the information and convert to a data frame row
    df_row = data.frame(
        var1 = result$var1,
        var2 = result$var2,
        p_value = result$test$p.value,
        odds_ratio = result$test$estimate, 
        total_expressed = result$total_expressed,
        degs = result$degs,
        kegg_true = result$kegg_TRUE,
        total_expressed_kegg = result$total_expressed_keg,
        fdr = p.adjust(result$test$p.value)
    )
    rownames(df_row) = NULL
    # Bind the data frame row to the final data frame
    final_df_greater_logfcgreater0 = rbind(final_df_greater_logfcgreater0, df_row)
}

final_df_greater_logfcgreater0 = final_df_greater_logfcgreater0[which(final_df_greater_logfcgreater0$fdr <= 0.05),]
final_df_greater_logfcgreater0$test = "greater_logfcgreater0"
final_df_greater_logfcgreater0$enrichment = "oxidative_phosphorylation"
saveRDS(final_df_greater_logfcgreater0, "oxidative_phosphorylation_result.RDS")
