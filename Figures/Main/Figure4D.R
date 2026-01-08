#################
library(dplyr)
library(ggplot2)
library(tidyr)

myres = readRDS("sen_sig_GO_by_features_7.8.24.RDS")

myres$significant = ifelse(myres$BH <= 0.05, TRUE, FALSE)
myres2 = myres[which(myres$BH <= 0.05),]
myres2 = as.data.frame(myres2)
# myres2$feature <- sub("_[^_]*$", "", myres2$feature)
myres2 = myres2[myres2$celltype != "noneu",]

myres2 = myres2[myres2$direction != "negative",]

res = data.frame(GO.ID = unique(myres2$GO.ID), mg= NA, exc= NA, oli = NA,int = NA)
for(x in unique(myres2$GO.ID)){
  sub = myres2[myres2$GO.ID == x,]
  tab= table(sub$celltype)
  res$mg[res$GO.ID == x] = tab["mg"]
  res$exc[res$GO.ID == x] = tab["exc"]
  res$oli[res$GO.ID == x] = tab["oli"]
  res$int[res$GO.ID == x] = tab["int"]

}

res$ncell_type = rowSums(res[,2:5]/res[,2:5], na.rm = TRUE)
res$mean_go_term = rowMeans(res[,2:5],na.rm = TRUE)

res = res[order(res$mean_go_term),]
res = res[order(res$ncell_type),]

###panel 1 One across celltypes bottom 30
myres2_bottom30 <- tail(res, 16)


myres_final = myres2[myres2$GO.ID %in% unique(myres2_bottom30$GO.ID), ]

myres_final$celltype = gsub("mg", "Microglia", myres_final$celltype)
myres_final$celltype = gsub("exc", "Excitatory Neurons", myres_final$celltype)
myres_final$celltype = gsub("oli", "Oligodendrocytes", myres_final$celltype)
myres_final$celltype = gsub("int", "Inhibitory Neurons", myres_final$celltype)
myres_final$celltype = gsub("opc", "Oligodendrocyte Progenitor Cells", myres_final$celltype)
myres_final$celltype = gsub("ast", "Astrocytes", myres_final$celltype)


myres_final = order[myres_final$fold_enrichment]


p <- ggplot(
  myres_final,
  aes(
    x = celltype,
    y = reorder(Term, -log10(BH), FUN = max),  # order by max significance
    fill = -log10(BH)
  )
) +
  geom_tile(color = "white") +  # white grid lines
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1,
    name = "-log10(FDR)"
  ) +
  labs(x = "Cell Type", y = "GO Term") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.background = element_rect(
      colour = "black",
      size = 0.5,
      fill = "transparent"
    )
  )

ggsave("sen_p_value_GOA.pdf", p, height = 10, width = 10)
