# Load necessary libraries
library(ggseg)
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape)
library(stringr)
library(tidyr) 

myres = readRDS("cor_smri_aucell_liv.RDS")
pm = readRDS("cor_upto5_smri_aucell.RDS")
par3 = readRDS("cor_5andup_smri_aucell.RDS")

myres$ds = "liv"
pm$ds = "upto5"
par3$ds = "fiveandup"


pm$assay = gsub("micro", "mg", pm$assay)
pm$assay = gsub("Astro", "ast", pm$assay)
pm$assay = gsub("OPC", "opc", pm$assay)
pm$assay = gsub("Oligo", "oli", pm$assay)
pm$assay = gsub("vas", "noneu", pm$assay)

par3$assay = gsub("micro", "mg", par3$assay)
par3$assay = gsub("Astro", "ast", par3$assay)
par3$assay = gsub("OPC", "opc", par3$assay)
par3$assay = gsub("Oligo", "oli", par3$assay)
par3$assay = gsub("vas", "noneu", par3$assay)

myres5 = rbind(myres, pm, par3)

myres5$assay = gsub("oli", "Oli", myres5$assay)
myres5$assay = gsub("exc", "Exc", myres5$assay)
myres5$assay = gsub("noneu", "NonNeu", myres5$assay)
myres5$assay = gsub("mg", "MG", myres5$assay)
myres5$assay = gsub("ast", "Ast", myres5$assay)
myres5$assay = gsub("opc", "OPC", myres5$assay)
myres5$assay = gsub("int", "Int", myres5$assay)

myres5.2 = myres5[which(myres5$pi1_sen_de > 0.01 & myres5$pi1_sc_de > 0.01),]
myres5.2 <- myres5.2 %>%
  mutate(new_column = str_extract(feature, "(?<=_)[^_]+$")) %>%
  mutate(new_column = ifelse(new_column %in% c("thickness", "area"), new_column, "volume"))

myres5.3 = myres5.2[which(myres5.2$p.value <= 0.05),]

myres_final = myres5.3[,c("feature","assay","spearmans", "ds")]


pi1_single = myres_final

# all = merge(myres2, pi1_single, by = c("feature","assay"),all = TRUE)

# pi1_single = myres_final
# pi1_single$top_score_sc = pi1_single$pi1Limma > 0

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

summary3$label = gsub("scale", "", summary3$label)

# Prepare the aseg atlas data
aseg2 <- as_ggseg_atlas(aseg) %>%
  na.omit() %>%
  distinct(label, .keep_all = TRUE)

# Merge the summary data with the atlas data
all2 <- merge(summary3, aseg2, by = "label")

# Group the data
all2_grouped <- all2 %>%
  group_by(assay, new_column)


all2_grouped2 = all2_grouped[all2_grouped$assay %in% c("MG", "Exc"),]


all2_grouped2$new_column = paste0(all2_grouped2$ds, all2_grouped2$assay) 
# Create the plot


b <- ggplot() +
  geom_brain(
    atlas = aseg,
    data  = all2_grouped2,
    aes(fill = spearmans)
  ) +
  facet_wrap(~new_column, ncol = 1) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-1, 1),                 # <- enforce range
    breaks = seq(-1, 1, by = 0.5),     # tidy legend ticks
    labels = scales::number_format(accuracy = 0.1),
    oob = scales::squish,              # clamp out-of-range values to [-1, 1]
    name = "Spearman's\n\u03C1",
    na.value = "gray80"
  ) +
  theme_void() +
  theme(
    legend.margin = margin(10, 10, 10, 10),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, hjust = 0.5),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt"),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10))

saveRDS(b, "volume_signal_pi1_ggseg_nd_b.rds")
