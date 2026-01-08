############
#volume
#############

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

# Filter and prepare the initial data
summary3 <- summary_d2f %>% 
  filter(new_column == "volume", assay != "NonNeu")

# Rename columns and clean up labels
summary3 <- summary3 %>%
  dplyr::rename(label = feature) %>%
  dplyr::mutate(label = gsub("_", "-", label))

# Prepare the aseg atlas data
aseg2 <- as_ggseg_atlas(aseg) %>%
  na.omit() %>%
  distinct(label, .keep_all = TRUE)

# Merge the summary data with the atlas data
all2 <- inner_join(summary3, aseg2, by = "label")

# Group the data
all2_grouped <- all2 %>%
  group_by(assay, new_column)

all2_grouped$new_column = gsub("volume", "Volume", all2_grouped$new_column)

g <- ggplot() +
  geom_brain(atlas = aseg, 
             data = all2_grouped, 
             aes(fill = pi1Limma)) +
  facet_wrap(~ assay, ncol = 1) +  # Facet by assay only
  scale_fill_gradient(low = "white", high = "blue", 
                      name = expression(pi[1]),  # Use expression for π₁ symbol
                      limits = c(0, 0.4), 
                      na.value = "white") +
  theme_void() +   
  theme(
    legend.margin = margin(10, 10, 10, 10),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16, hjust = 0.2),  # Center legend title
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt"),
    strip.text = element_text(size = 16, face = "bold")  # Increase facet label size
  ) +
  guides(
    fill = guide_colorbar(barwidth = 2, barheight = 10)  # Increase color bar size
  )

ggsave("signal_volume_pi1_ggseg.pdf",g) 
