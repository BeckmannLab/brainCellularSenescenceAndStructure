#################
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

myres = readRDS("sen_sig_GO_by_features_7.8.24.RDS")

myres$significant = ifelse(myres$BH <= 0.05, TRUE, FALSE)
myres2 = myres[which(myres$BH <= 0.05),]
myres2 = as.data.frame(myres2)
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


res2 = merge(unique(myres2[,c("GO.ID", "Term")]), res)

res2 = res2[order(res2$mean_go_term),]
res2 = res2[order(res2$ncell_type),]

###panel 1 One across celltypes bottom 30
myres2_bottom30 <- tail(res2, 30)


df_long <- myres2_bottom30 %>%
  pivot_longer(cols = c(mg, exc, oli, int), 
               names_to = "CellType", 
               values_to = "Count")



df_long <- df_long %>%
  group_by(Term) %>%
  mutate(Total = sum(Count)) %>%
  ungroup()


df_long$CellType = gsub("mg", "Microglia", df_long$CellType)
df_long$CellType = gsub("exc", "Excitatory Neurons", df_long$CellType)
df_long$CellType = gsub("oli", "Oligodendrocytes", df_long$CellType)
df_long$CellType = gsub("int", "Inhibitory Neurons", df_long$CellType)
df_long$CellType = gsub("opc", "Oligodendrocyte Progenitor Cells", df_long$CellType)
df_long$CellType = gsub("ast", "Astrocytes", df_long$CellType)


color_scheme <- c(
  "Microglia" = "#1b9e77",
  "Excitatory Neurons" = "#e7298a",
  "Oligodendrocytes" = "#66a61e",
  "Inhibitory Neurons" = "#d95f02",
  "Astrocytes" = "#7570b3",
  "Oligodendrocyte Progenitor Cells" = "#e6ab02"
)


df_long$Term = gsub("cell-cell adhesion via plasma-membrane a...","cell-cell adhesion via plasma-membrane adhesion molecules", df_long$Term)


d <- ggplot(df_long, aes(x = Count, y = reorder(Term, ncell_type), fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Count", y = "GO Term", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 12, r = 12, b = 12, l = 12),
    legend.box.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
  ) +
  scale_fill_manual(values = color_scheme) +
  guides(fill = guide_legend(nrow = 2))
  
ggsave("sen_30_GOA.pdf", d, height = 10, width = 15)

##panel 2 other in only one cell type. top 30

result <- do.call(rbind, 
                  lapply(split(myres2, myres2$celltype), 
                         function(x) x[order(x$BH, decreasing=FALSE)[1:5],]))

myres = readRDS("sen_sig_GO_by_features_7.8.24.RDS")


myres$significant = ifelse(myres$BH <= 0.05, TRUE, FALSE)
myres2 = myres[which(myres$BH <= 0.05),]
myres2 = as.data.frame(myres2)
myres2 = myres2[myres2$celltype != "noneu",]

myres2 = myres2[myres2$direction != "negative",]

res = data.frame(GO.ID = unique(myres2$GO.ID), mg= NA, exc= NA, oli = NA,int = NA, opc = NA)
for(x in unique(myres2$GO.ID)){
  sub = myres2[myres2$GO.ID == x,]
  tab= table(sub$celltype)
  res$mg[res$GO.ID == x] = tab["mg"]
  res$exc[res$GO.ID == x] = tab["exc"]
  res$oli[res$GO.ID == x] = tab["oli"]
  res$int[res$GO.ID == x] = tab["int"]
  res$opc[res$GO.ID == x] = tab["opc"]

}
res$ncell_type = rowSums(res[,2:5]/res[,2:5], na.rm = TRUE)
res$mean_go_term = rowMeans(res[,2:5],na.rm = TRUE)


res2 = res[which(res$ncell_type==1),"GO.ID"]

res2 = as.data.frame(res2)
all = merge(myres2, res2, by.x = "GO.ID", by.y = "res2")

result <- do.call(rbind, 
                  lapply(split(all, all$celltype), 
                         function(x) x[order(x$BH, decreasing=FALSE)[1:10],]))


result$Term = factor(result$Term, levels = result$Term[order(result$fold_enrichment, decreasing = F)] )

result$celltype = gsub("exc", "Exc", result$celltype)
result$celltype = gsub("int", "Int", result$celltype)
result$celltype = gsub("oli", "Oli", result$celltype)
result$celltype = gsub("int", "Int", result$celltype)
result$celltype = gsub("ast", "Ast", result$celltype)
result$celltype = gsub("opc", "OPC", result$celltype)
result$celltype = gsub("mg", "MG", result$celltype)


color_scheme <- c(
  "MG" = "#1b9e77",
  "Exc" = "#e7298a",
  "Oli" = "#66a61e",
  "Int" = "#d95f02",
  "Ast" = "#7570b3",
  "OPC" = "#e6ab02",
)


g = ggplot(result, aes(x = fold_enrichment, y = reorder(Term, fold_enrichment), fill = celltype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",  # This removes the legend
    plot.margin = margin(10, 10, 10, 10)  # Adjusted margins
  ) +
  scale_fill_manual(values = color_scheme) +
  labs(x = "Fold Enrichment", y = "GO Term") +  # Removed 'fill' label
  coord_cartesian(clip = "off")  # Prevent clipping of labels

ggsave("sen_30_GO.pdf", g, height = 10, width = 15)
