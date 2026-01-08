 
ml R/4.4.1
R

library(foreach)
library(dreamlet)
library(dplyr)
library(tidyr)
library(doParallel)
registerDoParallel(cores=10)


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
), V4 = c("Micro_TRUE", "L2-3_CUX2_TRUE", "L4_RORB_TRUE", "L5-6_THEMIS_TRUE",
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

##up to 5
resfinal=foreach(i = 1:nrow(var)) %do% {
      if(length(which(paste0("correct_edger_upto5_",var[i,6],"_",var[i,2], "_fit.RDS")==list.files("path/to/results",recursive=TRUE)))==0){
    pb = readRDS(paste0(folder,var[i,1])) #pulling results made in Herring_DE_setup.R code
    subset_age <- pb[, c("RL2103_ga22_v3", "RL2107_ga24_v3", "RL2121_ga34_v3", "RL1777_2d_v3", "RL1612_34d_v2", "RL2100_86d_v3", "RL2104_118d_v3", "RL2108_179d_v3", "RL2122_301d_v3", "RL2125_422d_v3","RL2105_627d_v3", "RL1786_2yr_v3", "RL1613_2yr_v2", "RL2129_3yr_v3", "RL2109_4yr_v3")]
    ncells = as.data.frame(int_colData(subset_age)$n_cells)
    ncells$celltype = rownames(ncells) 
    ncells2 = ncells %>%
      pivot_longer(!celltype,names_to = "id", values_to = "count")
    ncells2 = as.data.frame(ncells2)
    ncells2$id = gsub("\\.", "_", ncells2$id)
    colnames(ncells2) = c("sen_auc", "Sample_ID", "ncell")

    ##add ncells
    metadata(subset_age)$aggr_means$Indv_ID = metadata(subset_age)$aggr_means$Sample_ID 

    #define pairs
    ct.pairs <- c(var[i,4], var[i,5])

    ##run function
    try(fit <- dreamletCompareClusters_edgeR(subset_age, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    #pull
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    #save
    saveRDS(fit_top,paste0(var[i,6],"_",var[i,2], "_fit.RDS"))

}else{
    print(paste(i, "done"))
}}


###5andup
resfinal=foreach(i = 1:nrow(var)) %do% {
      if(length(which(paste0("correct_edger_5andup_",var[i,6],"_",var[i,2], "_fit.RDS")==list.files("/path/to/results/",recursive=TRUE)))==0){
    pb = readRDS(paste0(folder,var[i,1]))
    subset_age <- pb[, c( "RL2106_6yr_v3", "RL1614_8yr_v2", "RL2126_10yr_v3", "RL2110_10yr_v3","RL2127_12yr_v3", "RL2102_16yr_v3", "RL2123_20yr_v3", "RL2124_40yr_v3", "RL2128_20yr_v3", "RL2130_14yr_v3", "RL2131_17yr_v3", "RL2132_25yr_v3")]
    ncells = as.data.frame(int_colData(subset_age)$n_cells)
    ncells$celltype = rownames(ncells) 

    ncells2 = ncells %>%
      pivot_longer(!celltype,names_to = "id", values_to = "count")

    ncells2 = as.data.frame(ncells2)

    ncells2$id = gsub("\\.", "_", ncells2$id)

    colnames(ncells2) = c("sen_auc", "Sample_ID", "ncell")

    ##add ncells
    metadata(subset_age)$aggr_means$Indv_ID = metadata(subset_age)$aggr_means$Sample_ID 

  	#pairs
    ct.pairs <- c(var[i,4], var[i,5])

    ##run function
    try(fit <- dreamletCompareClusters_edgeR(subset_age, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
    #pull
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    #save
    saveRDS(fit_top,paste0(var[i,6],"_",var[i,2], "_fit.RDS"))
}else{
    print(paste(i, "done"))
}}


##################
#enrichmenet
#################

##upto5
all2 = c()
for (i in 1:nrow(var)) {
    file_path = paste0(var[i,6], "_", var[i,2], "_fit.RDS")  
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        df$celltype = var[i,6]
        df$test = var[i,2]
        df$symbol = rownames(df)
        all2 = rbind(df, all2)
    }
}

all2$sig = all2$FDR <= 0.05

#senescence genes
sen_genes <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", 
"CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", 
"GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", 
"LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", 
"MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", 
"STING1", "TGFB1", "TIMP2", "TNF", "TP53","ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
"BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", 
"CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", 
"CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", 
"CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", 
"CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", 
"FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", 
"HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", 
"IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", 
"IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", 
"INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", 
"MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", 
"NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", 
"PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", 
"SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", 
"TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", 
"TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")

sen_genes = as.data.frame(unique(sen_genes))
colnames(sen_genes) = "symbol"
dim(sen_genes)
# [1] 146   1

sen_genes$sen_gene = TRUE
intersect(all2$symbol, sen_genes$symbol)
sc_with_sen=merge(all2,sen_genes, by = "symbol", all.x = T)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)]= FALSE

var_e = structure(list(V2 = c("1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "1percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "5percent", "10percent", 
"10percent", "10percent", "10percent", "10percent", "10percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"20percent", "20percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "30percent", "30percent", "1percent", 
"1percent", "5percent", "5percent", "10percent", "10percent", 
"20percent", "20percent", "30percent", "30percent", "1percent", 
"10percent", "20percent", "1percent", "5percent", "5percent", 
"5percent", "10percent", "10percent", "10percent", "20percent", 
"20percent", "20percent", "30percent", "30percent", "30percent"
), V6 = c("micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", 
"LAMP5", "PN_dev", "micro", "L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", 
"L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2", "L4_RORB", 
"L5-6_THEMIS", "L5-6_TLE4", "PN_dev", "micro", "L2-3_CUX2", "L4_RORB", 
"L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "micro", "L2-3_CUX2", 
"L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "LAMP5", "PN_dev", "int", 
"exc", "int", "exc", "int", "exc", "int", "exc", "int", "exc", 
"vas", "vas", "vas", "Astro", "Oligo", "OPC", "Astro", "Oligo", 
"OPC", "Astro", "Oligo", "OPC", "Astro", "Oligo", "OPC", "Astro"
)), row.names = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 
12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 21L, 22L, 23L, 24L, 25L, 
26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 
39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L, 48L, 49L, 53L, 54L, 55L, 
56L, 57L, 58L, 59L, 60L, 61L, 62L, 63L, 64L, 65L), class = "data.frame")

#run enrichment and pull results
results = list()
for (i in 1:nrow(var_e)) {
    sc_subset = sc_with_sen[which(sc_with_sen$celltype == var_e[i,2] & sc_with_sen$test == var_e[i,1]), ]
    sig_table = table(sc_subset$sig)
    
    if (length(sig_table) > 1 && all(c(TRUE, FALSE) %in% names(sig_table))) {
        test = try(fisher.test(table(sc_subset$sen_gene, sc_subset$sig)))
        results[[i]] = list(
            var1 = var_e[i,1], 
            var2 = var_e[i,2], 
            table = table(sc_subset$sen_gene, sc_subset$sig), 
            test = test, 
            total_expressed = length(sc_subset$celltype), 
            degs = length(which(sc_subset$FDR <= 0.05)), 
            kegg_TRUE = length(which(sc_subset$sen_gene == TRUE & sc_subset$sig == TRUE)), 
            total_expressed_kegg = length(which(sc_subset$sen_gene == TRUE))
        )
    }
}

final_df = data.frame(var1 = character(),
                      var2 = character(),
                      p_value = numeric(),
                      odds_ratio= numeric(),
                      total_expressed = numeric(),
                      degs = numeric(),
                      kegg_true = numeric(),
                      total_expressed_kegg = numeric()
                      )

# Iterate through the results list and extract the necessary information
for (i in 1:length(results)) {
    result = results[[i]]
    
    # Extract the information and convert to a data frame row
    df_row = data.frame(
        var1 = result$var1,
        var2 = result$var2,
        p_value = result$test$p.value,
        odds_ratio = result$test$estimate,
        degs = result$degs
    )
    
    # Bind the data frame row to the final data frame
    final_df = rbind(final_df, df_row)
}


final_df[which(final_df$var2 == "Astro"),] #
final_df[which(final_df$var2 == "exc"),]
final_df[which(final_df$var2 == "int"),]
final_df[which(final_df$var2 == "micro"),]
final_df[which(final_df$var2 == "Oligo"),]
final_df[which(final_df$var2 == "OPC"),]
final_df[which(final_df$var2 == "vas"),]

##5andup
all2 = c()
for (i in 1:nrow(var)) {
    file_path = paste0(var[i,6], "_", var[i,2], "_fit.RDS")  
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        df$celltype = var[i,6]
        df$test = var[i,2]
        df$symbol = rownames(df)
        all2 = rbind(df, all2)
    }
}

all2$sig = all2$FDR <= 0.05

#senescence genes
sen_genes <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", 
"CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", 
"GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", 
"LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", 
"MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", 
"STING1", "TGFB1", "TIMP2", "TNF", "TP53","ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
"BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", 
"CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", 
"CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", 
"CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", 
"CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", 
"FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", 
"HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", 
"IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", 
"IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", 
"INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", 
"MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", 
"NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", 
"PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", 
"SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", 
"TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", 
"TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")

sen_genes = as.data.frame(unique(sen_genes))
colnames(sen_genes) = "symbol"
dim(sen_genes)
# [1] 146   1

sen_genes$sen_gene = TRUE
intersect(all2$symbol, sen_genes$symbol)
sc_with_sen=merge(all2,sen_genes, by = "symbol", all.x = T)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)]= FALSE


#run enrichment and pull results
results = list()
for (i in 1:nrow(var_e)) {
    sc_subset = sc_with_sen[which(sc_with_sen$celltype == var_e[i,2] & sc_with_sen$test == var_e[i,1]), ]
    sig_table = table(sc_subset$sig)
    
    if (length(sig_table) > 1 && all(c(TRUE, FALSE) %in% names(sig_table))) {
        test = try(fisher.test(table(sc_subset$sen_gene, sc_subset$sig)))
        results[[i]] = list(
            var1 = var_e[i,1], 
            var2 = var_e[i,2], 
            table = table(sc_subset$sen_gene, sc_subset$sig), 
            test = test, 
            total_expressed = length(sc_subset$celltype), 
            degs = length(which(sc_subset$FDR <= 0.05)), 
            kegg_TRUE = length(which(sc_subset$sen_gene == TRUE & sc_subset$sig == TRUE)), 
            total_expressed_kegg = length(which(sc_subset$sen_gene == TRUE))
        )
    }
}

final_df = data.frame(var1 = character(),
                      var2 = character(),
                      p_value = numeric(),
                      odds_ratio= numeric(),
                      total_expressed = numeric(),
                      degs = numeric(),
                      kegg_true = numeric(),
                      total_expressed_kegg = numeric()
                      )

# Iterate through the results list and extract the necessary information
for (i in 1:length(results)) {
    result = results[[i]]
    
    # Extract the information and convert to a data frame row
    df_row = data.frame(
        var1 = result$var1,
        var2 = result$var2,
        p_value = result$test$p.value,
        odds_ratio = result$test$estimate,
        degs = result$degs
    )
    
    # Bind the data frame row to the final data frame
    final_df = rbind(final_df, df_row)
}


final_df[which(final_df$var2 == "Astro"),] #
final_df[which(final_df$var2 == "exc"),]
final_df[which(final_df$var2 == "int"),]
final_df[which(final_df$var2 == "micro"),]
final_df[which(final_df$var2 == "Oligo"),]
final_df[which(final_df$var2 == "OPC"),]
final_df[which(final_df$var2 == "vas"),]

