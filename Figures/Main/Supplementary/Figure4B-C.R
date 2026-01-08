ml R/4.4.1
R

library(SingleCellExperiment)
library(dreamlet)
library(muscat)
library(data.table)
library(assertthat)
library(foreach)
library(qvalue)

ast_fit = readRDS("ast_fit_ncells.RDS")
ast_point4 = readRDS("correct_edger_ast_10percent_fit.RDS")

data_ast <- topTable(ast_fit,coef = "compare", num=Inf)

cells = c("data_ast")


all2_voom = c()
for (x in cells){
  df = get(x)
  df$celltype = gsub("data_", "", x)
  all2_voom = rbind(df,all2_voom)
}

cells = c("ast_point4")

all2_edger = c()
for (x in cells){
  df = get(x)
  df$celltype = gsub("data_", "", x)
  all2_edger = rbind(df,all2_edger)
}


all2_voom$sig = all2_voom$adj.P.Val <= 0.05
all2_edger$sig = all2_edger$FDR <= 0.05

all3 = all2_voom[order(all2_voom$AveExpr, decreasing = F),]

#######
#limma
#########
all3$celltype = "dreamletCompareClusters"

y_limit <- max(abs(all3$logFC)) * 1.0 
limma_one <- ggplot(all3) +
  aes(x = AveExpr, y = logFC) + 
  geom_point(pch = 1, size = 3, alpha = 0.4, aes(col = sig), show.legend = TRUE) +
  scale_color_manual(
    name = expression("FDR" <= 0.05),
    values = c("#d80f8c", "#2cace2", "#999999")
  ) +
  geom_density2d(show.legend = FALSE, alpha = 0.4, colour = "black") +
  geom_smooth(show.legend = FALSE, colour = "black", se = FALSE) +
  scale_size(
    name = "FDR (-log10)",
    limits = c(0, 10), range = c(0, 0.5),
    guide = "none"  # Hides the size legend
  ) +
  coord_cartesian(ylim = c(-1, 1) * y_limit) +
  theme_minimal() + 
  facet_wrap(~celltype) +
  labs(
    x = "LogCPM",
    y = "LogFC"
  ) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15),
    legend.position = "bottom",
    legend.background = element_rect(color = "black", fill = "white", size = 0.5),
    legend.box.background = element_rect(color = "black", size = 0.5)
  )

        # theme(legend.position = "none") 


#################
#edgeR
#######################
all2_edger$celltype = "dreamletCompareClusters with negative binomial"

 y_limit <- max(abs(all2_edger$logFC)) * 1.0 
egder_one = ggplot(all2_edger) +
  aes(x = logCPM, y = logFC) + 
  geom_point(pch = 1, size = 3, alpha = 0.4, aes(col = sig), show.legend = TRUE) +  # Show legend
  scale_color_manual(
    values = c("#d80f8c", "#2cace2", "#999999"),
    name = expression("FDR" <= 0.05)  # Add a legend title
  ) +
  geom_density2d(alpha = 0.4, colour = "black", show.legend = FALSE) +
  geom_smooth(colour = "black", se = FALSE, show.legend = FALSE) +
  scale_size(
    name = "FDR (-log10)",
    limits = c(0, 10), range = c(0, 0.5), guide = "none"
  ) +
  coord_cartesian(ylim = c(-1, 1) * y_limit) +
  theme_minimal() + 
  facet_wrap(~celltype) +
  labs(
    x = "LogCPM",
    y = "LogFC"
  ) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15),
    legend.position = "bottom",
    legend.background = element_rect(color = "black", fill = "white", size = 0.5),
    legend.box.background = element_rect(color = "black", size = 0.5)
  )

library(patchwork)
method_compare = limma_one  /egder_one + theme(legend.position = "bottom") 

ggsave("edgeR_vs_limma.pdf", method_compare, width = 10, height = 10)
