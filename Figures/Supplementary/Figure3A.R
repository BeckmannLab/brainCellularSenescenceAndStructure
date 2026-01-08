#################
library(dplyr)
library(ggplot2)
library(tidyr)

myres = readRDS("new_annotation_bulk_pi1_sig_GO_by_features.RDS")
myres$feature <- sub("_[^_]*$", "", myres$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")



# Remove rows with specified row names
myres <- myres[!myres$feature %in% rows_to_remove, ]

myres$significant = ifelse(myres$BH <= 0.05, TRUE, FALSE)
myres2 = myres[which(myres$BH <= 0.05),]
myres2 = as.data.frame(myres2)
myres2$feature <- sub("_[^_]*$", "", myres2$feature)
myres2 = myres2[myres2$celltype != "NonNeu",]

res = data.frame(GO.ID = unique(myres2$GO.ID), mg= NA, exc= NA, oli = NA, OPC = NA, int = NA, ast = NA)
for(x in unique(myres2$GO.ID)){
  sub = myres2[myres2$GO.ID == x,]
  tab= table(sub$celltype)
  res$mg[res$GO.ID == x] = tab["MG"]
  res$exc[res$GO.ID == x] = tab["Exc"]
  res$oli[res$GO.ID == x] = tab["Oli"]
  res$OPC[res$GO.ID == x] = tab["OPC"]
  res$int[res$GO.ID == x] = tab["Int"]
  res$ast[res$GO.ID == x] = tab["Ast"]
}

res$ncell_type = rowSums(res[,2:7]/res[,2:7], na.rm = TRUE)
res$mean_go_term = rowMeans(res[,2:7],na.rm = TRUE)


res2 = merge(unique(myres2[,c("GO.ID", "Term")]), res)
res2 = res2[order(res2$mean_go_term),]
res2 = res2[order(res2$ncell_type),]

myres2_bottom30 <- tail(res2, 30)

df_long <- myres2_bottom30 %>%
  pivot_longer(cols = c(mg, exc, oli, OPC, int, ast), 
               names_to = "CellType", 
               values_to = "Count")

# Calculate total count for each term
df_long <- df_long %>%
  group_by(Term) %>%
  mutate(Total = sum(Count)) %>%
  ungroup()

df_long$CellType2 = df_long$CellType
df_long$CellType2 = gsub("mg", "Microglia", df_long$CellType2)
df_long$CellType2 = gsub("exc", "Excitatory Neurons", df_long$CellType2)
df_long$CellType2 = gsub("oli", "Oligodendrocytes", df_long$CellType2)
df_long$CellType2 = gsub("int", "Inhibitory Neurons", df_long$CellType2)
df_long$CellType2 = gsub("ast", "Astrocytes", df_long$CellType2)
df_long$CellType2 = gsub("OPC", "Oligodendrocyte progenitor cells", df_long$CellType2)


color_scheme <- c(
  "Microglia" = "#1b9e77",
  "Excitatory Neurons" = "#e7298a",
  "Oligodendrocytes" = "#66a61e",
  "Inhibitory Neurons" = "#d95f02",
  "Astrocytes" = "#7570b3",
  "Oligodendrocyte progenitor cells" = "#e6ab02"
)

# Create the stacked bar plot
g = ggplot(df_long, aes(x = Count, y = reorder(Term, Total), fill = CellType2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.5) +  # Increased size
  labs(x = "Count", y = "GO Term", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 16),       # Slightly bigger
    axis.text.x = element_text(size = 16),       # Slightly bigger
    axis.title.y = element_text(size = 18),      # Slightly bigger
    axis.title.x = element_text(size = 18),      # Slightly bigger
    legend.position = c(0.95, 0.1),
    legend.justification = c(1, 0),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size = 14),       # Slightly bigger
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 16),      # Slightly bigger
    legend.title.align = 0.5
  ) +
  scale_fill_manual(values = color_scheme) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, max(df_long$Total), by = 20)
  ) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(ncol = 1))  # One column legend

ggsave("single_GO_top30.pdf", g, height = 10, width = 15,limitsize = FALSE)
