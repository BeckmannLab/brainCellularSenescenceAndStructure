#############
# figure 4B
#############

library(ggplot2)
library(data.table)
library(ggrepel)
library(dplyr)

var  <- read.delim(
  "var_de_aucell.txt",
  header = FALSE
)
var2 <- read.delim(
  "var_de_aucell_10.txt",
  header = FALSE
)

var <- rbind(var, var2)

all2 <- c()
for (i in 1:nrow(var)) {
  df <- readRDS(paste0(
    "correct_edger_",
    var[i, 6], "_", var[i, 2], "_fit.RDS"
  ))
  df$celltype <- var[i, 6]
  df$test     <- var[i, 2]
  df$symbol   <- rownames(df)
  all2        <- rbind(df, all2)
}

liv_aucell <- all2
liv_aucell$celltype <- gsub("noneu", "nonneu", liv_aucell$celltype)

all2$sig           <- all2$FDR <= 0.05
all2               <- as.data.table(all2)

all3 <- all2[
  (celltype == "oli" & test == "10percent") |
  (celltype == "exc" & test == "5percent")  |
  (celltype == "int" & test == "5percent")  |
  (celltype == "mg"  & test == "30percent") |
  (celltype == "ast" & test == "1percent")  |
  (celltype == "opc" & test == "30percent")
]

all3$logFC_sign <- factor(
  ifelse(all3$logFC > 0, "Positive", "Negative"),
  levels = c("Positive", "Negative")
)

all3$new_celltype <- all3$celltype
all3$new_celltype <- gsub("ast", "Astrocytes",                      all3$new_celltype)
all3$new_celltype <- gsub("exc", "Excitatory Neurons",              all3$new_celltype)
all3$new_celltype <- gsub("int", "Inhibitory Neurons",              all3$new_celltype)
all3$new_celltype <- gsub("mg",  "Microglia",                       all3$new_celltype)
all3$new_celltype <- gsub("oli", "Oligodendrocytes",                all3$new_celltype)
all3$new_celltype <- gsub("opc", "Oligodendrocyte progenitor cells", all3$new_celltype)

b <- ggplot(all3) +
  aes(x = logFC, y = -log10(PValue)) +
  geom_point(aes(color = case_when(
    !sig          ~ "Not Significant",
    sig & logFC > 0 ~ "Upregulated",
    sig & logFC < 0 ~ "Downregulated"
  )), shape = 16) +
  facet_wrap(~new_celltype, scales = "free") +
  scale_color_manual(
    values = c(
      "Not Significant" = "grey",
      "Upregulated"     = "red",
      "Downregulated"   = "blue"
    ),
    name = "Significance"
  ) +
  theme_bw() +
  labs(
    x     = expression("Log"[2] * " Fold Change"),
    y     = expression("-" * "Log"[10] * "(P-value)"),
    color = "Significant\n(FDR <= 0.05)"
  ) +
  theme(
    axis.text.x       = element_text(size = 20),
    axis.text.y       = element_text(size = 20),
    axis.title.x      = element_text(size = 24),
    axis.title.y      = element_text(size = 24),
    legend.text       = element_text(size = 22),
    legend.title      = element_text(size = 24),
    strip.text        = element_text(size = 20),
    legend.key.size   = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(1, "cm"),
    legend.position   = "bottom",
    legend.background = element_rect(color = "black", fill = "white")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))


ggsave(
  "volcano_plot_sen_DE.pdf",
  b, width = 12, height = 11
)
