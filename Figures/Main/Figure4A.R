library(seriation)
library(ggplot2)
library(data.table)
library(dplyr)
library(qvalue)
library(limma)

genes_42 <- c(
  "4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", 
  "CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", 
  "GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", 
  "LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", 
  "MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", 
  "STING1", "TGFB1", "TIMP2", "TNF", "TP53"
)

genes_125 <- c(
  "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", 
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
  "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2"
)

all <- c(genes_42, genes_125)
all <- unique(all)

final_df <- c()
new_df <- c()

var <- read.delim("var_de_aucell.txt", header = FALSE)

V1 <- rep("lbp_06.21.24_pseudobulk_10percent.RDS", 7)
V2 <- rep("10percent", 7)
V3 <- rep("pb_10percent", 7)

V4 <- c("Ast_TRUE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE")
V5 <- c("Ast_FALSE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE")
V6 <- c("ast", "exc", "int", "mg", "noneu", "oli", "opc")
V7 <- rep("0.4", 7)

var2 <- as.data.frame(cbind(V1, V2, V3, V4, V5, V6, V7))
var3 <- rbind(var, var2)
var <- var3

for (i in 1:nrow(var)) {
  df <- readRDS(
    paste0(
      "correct_edger_",
      var[i, 6], "_", var[i, 2], "_fit.RDS"
    )
  )
  
  new_df$sig        <- length(which(df$FDR <= 0.05))
  new_df$sen_genes  <- length(intersect(rownames(df), all))
  sig               <- length(which(df$FDR <= 0.05))
  new_df$total      <- nrow(df)
  total             <- nrow(df)
  new_df$percentage <- sig / total
  new_df$celltype   <- var[i, 6]
  new_df$test       <- var[i, 2]
  
  length_up   <- sum(df$logFC > 0)
  length_down <- sum(df$logFC < 0)
  
  new_df$percentage_up   <- length_up / total
  new_df$percentage_down <- length_down / total
  new_df$pi1             <- 1 - propTrueNull(df$PValue)
  
  new_df   <- as.data.frame(new_df)
  final_df <- rbind(final_df, new_df)
}

final_df3 <- final_df

final_df <- readRDS(
  "clean_across_percent_fishers_results_9_21_24.RDS"
)

colnames(final_df)[1] <- "test"
colnames(final_df)[2] <- "celltype"

all_percentages_celltype <- merge(final_df, final_df3, by = c("test", "celltype"))

all_percentages_celltype$log_oddsratio <- log2(all_percentages_celltype$odds_ratio)

max_non_inf <- max(
  all_percentages_celltype$log_oddsratio[is.finite(all_percentages_celltype$log_oddsratio)], 
  na.rm = TRUE
)

all_percentages_celltype$log_oddsratio[is.infinite(all_percentages_celltype$log_oddsratio)] <- 
  max_non_inf + 1

all_percentages_celltype <- as.data.table(all_percentages_celltype)

seriation_result <- seriate(
  dist(all_percentages_celltype %>% select(celltype, test, log_oddsratio)), 
  method = "OLO"
)

ordered_indices <- get_order(seriation_result)
all_ordered     <- all_percentages_celltype[ordered_indices, ]

all_ordered$test <- factor(
  all_ordered$test, 
  levels = c("1percent", "5percent", "10percent", "20percent", "30percent")
)

all_ordered2 <- all_ordered[all_ordered$celltype != "noneu", ]

all_ordered2$new_celltype <- all_ordered2$celltype
all_ordered2$new_celltype <- gsub("ast", "Astrocytes", all_ordered2$new_celltype)
all_ordered2$new_celltype <- gsub("exc", "Excitatory Neurons", all_ordered2$new_celltype)
all_ordered2$new_celltype <- gsub("int", "Inhibitory Neurons", all_ordered2$new_celltype)
all_ordered2$new_celltype <- gsub("mg", "Microglia", all_ordered2$new_celltype)
all_ordered2$new_celltype <- gsub("oli", "Oligodendrocytes", all_ordered2$new_celltype)
all_ordered2$new_celltype <- gsub(
  "opc", "Oligodendrocyte\nprogenitor cells", 
  all_ordered2$new_celltype
)

all_ordered2$test <- gsub("percent", "%", all_ordered2$test)
all_ordered2$test <- factor(
  all_ordered2$test, 
  levels = c("1%", "5%", "10%", "20%", "30%")
)

upper_limit   <- ceiling(max(all_ordered2$log_oddsratio, na.rm = TRUE))
legend_breaks <- c(0, 5, 10, upper_limit)
legend_labels <- c("0", "5", "10", "Inf")

d <- ggplot(all_ordered2, aes(x = new_celltype, y = test, fill = log_oddsratio)) +
  geom_tile(col = "white") +
  geom_text(
    aes(
      label = sprintf(
        "%.1f\n(p=%s)", 
        odds_ratio, 
        formatC(p_value, format = "e", digits = 1)
      )
    ),
    size = 8, 
    lineheight = 0.8
  ) +
  scale_fill_gradient(
    low    = "white", 
    high   = "red",
    limits = c(0, upper_limit),
    breaks = legend_breaks,
    labels = legend_labels,
    name   = expression(atop("Log"["2"], "(Odds Ratio)"))
  ) +
  labs(
    x = "Cell Type",
    y = "Senescence cell proportion threshold"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 25),
    axis.title.x     = element_text(size = 25),
    axis.title.y     = element_text(size = 27),
    legend.title     = element_text(size = 22, margin = margin(b = 10), hjust = 0.5),
    legend.title.align = 0.5,
    legend.text      = element_text(size = 20),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.justification = "center",
    plot.title       = element_text(size = 24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      frame.colour  = "black",
      frame.linewidth = 0.2,
      ticks.colour  = "black",
      barheight     = 15,
      barwidth      = 1.5
    )
  )

ggsave(
  "DE_oddsratio_heatmap_sen.pdf", 
  d, width = 10, height = 10
)
