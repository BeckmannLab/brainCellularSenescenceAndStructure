library(seriation)
library(ggplot2)
library(data.table)
library(dplyr)


##upto5
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

all2 = c()
for (i in 1:nrow(var)) {
    file_path = paste0("/correct_edger_upto5_", var[i,6], "_", var[i,2], "_fit.RDS")  
    if (file.exists(file_path)) {
        df = readRDS(file_path)
        df$celltype = var[i,6]
        df$test = var[i,2]
        df$symbol = rownames(df)
        all2 = rbind(df, all2)
    }
}

all2$sig = all2$FDR <= 0.05


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

var = structure(list(V1 = c("lbp_06.21.24_pseudobulk_10percent.RDS", 
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
    0.4, 0.4, 0.4, 0.4, 0.4)), row.names = c(NA, 7L), class = "data.frame")

head(var)

results = list()
for (i in 1:nrow(var)) {
    sc_subset = sc_with_sen[which(sc_with_sen$celltype == var[i,2] & sc_with_sen$test == var[i,1]), ]
    sig_table = table(sc_subset$sig)
    
    if (length(sig_table) > 1 && all(c(TRUE, FALSE) %in% names(sig_table))) {
        test = try(fisher.test(table(sc_subset$sen_gene, sc_subset$sig)))
        results[[i]] = list(
            var1 = var[i,1], 
            var2 = var[i,2], 
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


final_df2= final_df

final_df3= final_df
final_df3$ds = "uptofive"


final_df2$test = gsub("1percent", "1", final_df2$test)
final_df2$test = gsub("10percent", "10", final_df2$test)
final_df2$test = gsub("20percent", "20", final_df2$test)
final_df2$test = gsub("30percent", "30", final_df2$test)
final_df2$test = gsub("5percent", "5", final_df2$test)

final_df2$test = factor(final_df2$test, levels = c("1", "5", "10", "20", "30"))


colnames(final_df)[1] = "test"
colnames(final_df)[2] = "celltype"

final_df$log_oddsratio = log2(final_df$odds_ratio)

max_non_inf <- max(final_df$log_oddsratio[is.finite(final_df$log_oddsratio)], na.rm = TRUE)
final_df$log_oddsratio[is.infinite(final_df$log_oddsratio)] <- max_non_inf + 1

final_df = as.data.table(final_df)


seriation_result <- seriate(dist(final_df %>% select(celltype, test, log_oddsratio)), method = "OLO")

# Get the ordering
ordered_indices <- get_order(seriation_result)

# Order the data frame
all_ordered <- final_df[ordered_indices, ]

all_ordered$test = factor(all_ordered$test, levels = c("1percent", "5percent", "10percent", "20percent", "30percent"))

all_ordered2 = all_ordered[which(all_ordered$celltype!= "vas"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "LAMP5"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "L2-3_CUX2"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "L4_RORB"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "L5-6_THEMIS"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "L5-6_TLE4"),]
all_ordered2 = all_ordered2[which(all_ordered2$celltype!= "PN_dev"),]


all_ordered2$celltype = gsub("micro", "Microglia", all_ordered2$celltype)

all_ordered2$celltype = gsub("Astro", "Astrocytes", all_ordered2$celltype)
all_ordered2$celltype = gsub("exc", "Excitatory Neurons", all_ordered2$celltype)
all_ordered2$celltype = gsub("int", "Inhibitory Neurons", all_ordered2$celltype)

all_ordered2$celltype = gsub("Oligo", "Oligodendrocytes", all_ordered2$celltype)
all_ordered2$celltype = gsub("OPC", "Oligodendrocyte\nProgenitor Cells", all_ordered2$celltype)

all_ordered2$test = gsub("percent", "%",all_ordered2$test  )

all_ordered2$test = factor(all_ordered2$test, levels = c("1%", "5%", "10%", "20%", "30%"))
##plot
max_val <- max(all_ordered2$log_oddsratio, na.rm = TRUE)

upper_limit <- ceiling(max(all_ordered2$log_oddsratio, na.rm = TRUE))
legend_breaks <- c(0, 5, 10, upper_limit)
legend_labels <- c("0", "5", "10", "Inf")

##plot
d <- ggplot(all_ordered2, aes(x = celltype, y = test, fill = log_oddsratio)) +
  geom_tile(col = "white") +
  geom_text(aes(label = sprintf("%.1f\n(p=%s)", odds_ratio, 
                                formatC(p_value, format = "e", digits = 1))), 
            size = 8, lineheight = 0.8) +
  scale_fill_gradient2(
    low = "blue",     
    mid = "white",    
    high = "red",     
    midpoint = 0,     
    name = expression(atop("Log"["2"], "(Odds Ratio)"))
  ) +
  labs(
    x = "Cell Type",
    y = "Senescence cell proportion threshold"
  ) +
  theme_bw() +  
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 27),
    legend.title = element_text(size = 22, margin = margin(b = 10), hjust = 0.5),
    legend.title.align = 0.5,
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = "center",
    plot.title = element_text(size = 24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  coord_fixed() + 
  guides(fill = guide_colorbar(
    frame.colour = "black", 
    frame.linewidth = 0.2,
    ticks.colour = "black",
    barheight = 15,
    barwidth = 1.5
  ))


ggsave("ND_upto5_DE_oddsratio_heatmap_sen.pdf", d, width = 10, height = 10)
