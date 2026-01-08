# Load necessary libraries
library(ggseg)
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape)

pi1_single = readRDS("single_pi1_per_feature_celltype_09.16.24.RDS")

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
             "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
             "WM_hypointensities", "non_WM_hypointensities",
             "MaskVol", "BrainSegVol_to_eTIV",
             "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles", "SurfaceHoles")

pi1_single2 <- pi1_single[!pi1_single$feature %in% rows_to_remove, ]

summary_d2f <- pi1_single2 %>%
  mutate(
    new_column = case_when(
      grepl("thickness", feature) ~ "thickness",
      grepl("area", feature) ~ "area",
      TRUE ~ "volume"
    ),
    feature = gsub("thickness_|area_", "", feature)
  )


summary3 = summary_d2f[which(summary_d2f$new_column != "volume"),]
summary4 = summary3[which(summary3$assay != "NonNeu"),]
summary4$feature = gsub("_thickness|_area", "", summary4$feature)

colnames(summary4)[2] = "label"

summary_df3 = summary4[,c("label", "pi1Limma", "assay", "new_column")]

dk2 <- as_ggseg_atlas(dk)
dk2 <- na.omit(dk2)
dk2 <- dk2[!duplicated(dk2$label), ]

all2 = merge(summary_df3, dk2, by = "label")

# Prepare the summary data
summary_df4 <- summary_df3[!is.na(match(summary_df3$label, dk2$label)), ]

all2 <- merge(summary_df4, dk2, by = "label", all.x = TRUE)

all2_grouped <- all2 %>%
  group_by(assay, new_column)

all2_grouped$new_column = gsub("area", "Area", all2_grouped$new_column)
all2_grouped$new_column = gsub("thickness", "Thickness", all2_grouped$new_column)

all2_grouped$new_celltype = all2_grouped$assay
all2_grouped$new_celltype = gsub("MG", "Microglia", all2_grouped$new_celltype)
all2_grouped$new_celltype = gsub("Exc", "Excitatory Neurons", all2_grouped$new_celltype)
all2_grouped$new_celltype = gsub("Int", "Inhibitory Neurons", all2_grouped$new_celltype)
all2_grouped$new_celltype = gsub("Oli", "Oligodendrocytes", all2_grouped$new_celltype)
all2_grouped$new_celltype = gsub("Ast", "Astrocytes", all2_grouped$new_celltype)
all2_grouped$new_celltype = gsub("OPC", "Oligodendrocyte\nprogenitor cells", all2_grouped$new_celltype)

g =ggplot() +
  geom_brain(
    atlas = dk, 
    data = all2_grouped, 
    aes(fill = pi1Limma)
  ) +
  facet_grid(assay ~ new_column) +
  scale_fill_gradient(
    low = "white", 
    high = "blue", 
    name = "pi1", 
    na.value = "white",
    limits = c(0, 0.4),
    oob = scales::squish
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(10, 10, 10, 10)
  )

 
ggsave("signal_thickness_area_pi1_ggseg.pdf",g, width = 15, height = 10) 
