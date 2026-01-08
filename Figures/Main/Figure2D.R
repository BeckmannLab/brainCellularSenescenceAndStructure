#################
# Fig 2D
#################

library(ggplot2)
library(tidyr)

##-------------------------------
## Load and filter results
##-------------------------------

myres <- readRDS("new_annotation_bulk_pi1_sig_GO_by_features.RDS")

# Strip final suffix from feature name
myres$feature <- sub("_[^_]*$", "", myres$feature)

# Remove unwanted features
rows_to_remove <- c(
  "Left_WM_hypointensities", "Right_WM_hypointensities",
  "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
  "WM_hypointensities", "non_WM_hypointensities",
  "MaskVol", "BrainSegVol_to_eTIV",
  "MaskVol_to_eTIV", "SurfaceHoles",
  "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles"
)

myres <- myres[!(myres$feature %in% rows_to_remove), ]

# Mark significance and keep only BH <= 0.05 and non-NonNeu
myres$significant <- myres$BH <= 0.05

myres2 <- myres[myres$BH <= 0.05, ]
myres2 <- myres2[myres2$celltype != "NonNeu", ]

##-------------------------------
## Count significant GO terms per cell type
##-------------------------------

cell_types <- c("MG", "Exc", "Oli", "OPC", "Int", "Ast")

res <- data.frame(
  GO.ID = unique(myres2$GO.ID),
  mg  = NA_real_,
  exc = NA_real_,
  oli = NA_real_,
  OPC = NA_real_,
  int = NA_real_,
  ast = NA_real_
)

for (x in unique(myres2$GO.ID)) {
  sub <- myres2[myres2$GO.ID == x, ]
  tab <- table(sub$celltype)

  res$mg[res$GO.ID == x]  <- tab["MG"]
  res$exc[res$GO.ID == x] <- tab["Exc"]
  res$oli[res$GO.ID == x] <- tab["Oli"]
  res$OPC[res$GO.ID == x] <- tab["OPC"]
  res$int[res$GO.ID == x] <- tab["Int"]
  res$ast[res$GO.ID == x] <- tab["Ast"]
}

# Number of cell types where GO term is present (same logic as original x/x trick)
res$ncell_type   <- rowSums(res[, 2:7] / res[, 2:7], na.rm = TRUE)
res$mean_go_term <- rowMeans(res[, 2:7], na.rm = TRUE)

##-------------------------------
## Attach GO term names and pick top 10
##-------------------------------

res2 <- merge(
  unique(myres2[, c("GO.ID", "Term")]),
  res,
  by = "GO.ID"
)

# Same two-step ordering as in original code
res2 <- res2[order(res2$mean_go_term), ]
res2 <- res2[order(res2$ncell_type), ]

# Take last 10 (those with highest ncell_type after ordering)
myres2_bottom30 <- tail(res2, 10)

##-------------------------------
## Long format + totals per term
##-------------------------------

df_long <- pivot_longer(
  myres2_bottom30,
  cols      = c(mg, exc, oli, OPC, int, ast),
  names_to  = "CellType",
  values_to = "Count"
)

# Total count per GO term (same behaviour as sum(Count) with default NA handling)
total_per_term <- tapply(df_long$Count, df_long$Term, sum)
df_long$Total  <- total_per_term[df_long$Term]

# Pretty cell type labels
df_long$CellType2 <- df_long$CellType
df_long$CellType2 <- gsub("mg",  "Microglia",                       df_long$CellType2)
df_long$CellType2 <- gsub("exc", "Excitatory Neurons",              df_long$CellType2)
df_long$CellType2 <- gsub("oli", "Oligodendrocytes",                df_long$CellType2)
df_long$CellType2 <- gsub("int", "Inhibitory Neurons",              df_long$CellType2)
df_long$CellType2 <- gsub("ast", "Astrocytes",                      df_long$CellType2)
df_long$CellType2 <- gsub("OPC", "Oligodendrocyte progenitor cells", df_long$CellType2)

##-------------------------------
## Color scheme
##-------------------------------

color_scheme <- c(
  "Microglia"                       = "#1b9e77",
  "Excitatory Neurons"             = "#e7298a",
  "Oligodendrocytes"               = "#66a61e",
  "Inhibitory Neurons"             = "#d95f02",
  "Astrocytes"                     = "#7570b3",
  "Oligodendrocyte progenitor cells" = "#e6ab02"
)

##-------------------------------
## Plot
##-------------------------------

g <- ggplot(df_long, aes(x = Count, y = reorder(Term, Total), fill = CellType2)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 8
  ) +
  labs(
    x    = "Count",
    y    = "GO Term",
    fill = "Cell Type"
  ) +
  scale_fill_manual(values = color_scheme) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, max(df_long$Total, na.rm = TRUE), by = 20)
  ) +
  theme_minimal() +
  theme(
    axis.text        = element_text(size = 25),
    axis.title       = element_text(size = 28),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.margin    = margin(6, 6, 6, 6),
    legend.text      = element_text(size = 25),
    legend.title     = element_text(size = 27),
    legend.title.align = 0.5,
    legend.background  = element_rect(fill = "white"),
    legend.key.width   = unit(1.5, "cm"),
    legend.key.height  = unit(0.8, "cm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 6,
      override.aes = list(size = 5)
    )
  )

ggsave("single_GO_top10.pdf",g, height = 7, width = 20, limitsize = FALSE)
