##########
# Fig 2C
##########

# setup 
library(dreamlet)
library(data.table)
library(qvalue)
library(ggplot2)
library(ggthemes)
library(forcats)
library(scales)

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

pi1_single <- readRDS("single_pi1_per_feature_celltype_09.16.24.RDS")

##------------------------------------------------------------
## Filter features and assays
##------------------------------------------------------------

rows_to_remove <- c(
  "Left_WM_hypointensities", "Right_WM_hypointensities",
  "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
  "WM_hypointensities", "non_WM_hypointensities",
  "MaskVol", "BrainSegVol_to_eTIV",
  "MaskVol_to_eTIV", "SurfaceHoles",
  "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles"
)

pi1_single <- pi1_single[!(pi1_single$feature %in% rows_to_remove), ]
pi1_single <- pi1_single[pi1_single$assay != "NonNeu", ]

pi1_single2 <- pi1_single

##------------------------------------------------------------
## Add type (area / thickness / volume) – same logic as case_when
##------------------------------------------------------------

pi1_single2$type <- "volume"
idx_area      <- grepl("area",      pi1_single2$feature, ignore.case = TRUE)
idx_thickness <- grepl("thickness", pi1_single2$feature, ignore.case = TRUE)

pi1_single2$type[idx_area]      <- "area"
pi1_single2$type[idx_thickness] <- "thickness"

##------------------------------------------------------------
## Define PFC regions and flag PFC features
## (same grep logic as original)
##------------------------------------------------------------

pfc <- data.frame(
  feature  = c(
    "caudalmiddlefrontal", "lateralorbitofrontal", "medialorbitofrontal",
    "parsopercularis", "parsorbitalis", "parstriangularis",
    "rostralmiddlefrontal", "superiorfrontal", "frontalpole"
  ),
  feature1 = NA_character_,
  feature2 = NA_character_
)

for (feature in pfc$feature) {
  tmp <- grep(feature, pi1_single2$feature, value = TRUE)
  # same assumption as original: at least 2 matches per pattern
  pfc$feature1[pfc$feature == feature] <- tmp[1]
  pfc$feature2[pfc$feature == feature] <- tmp[2]
}

pfc2 <- data.frame(
  feature = c(pfc$feature1, pfc$feature2),
  inPFC   = TRUE,
  stringsAsFactors = FALSE
)

pfc_single <- merge(pfc2, pi1_single2, by = "feature", all.y = TRUE)
pfc_single$inPFC[is.na(pfc_single$inPFC)] <- FALSE

##------------------------------------------------------------
## Proportion of signatures with pi1 > 0.01 per assay
## (does not affect p-values, but kept for completeness)
##------------------------------------------------------------

res <- data.frame(
  assay      = unique(pi1_single2$assay),
  proportion = NA_real_
)

for (a in unique(pi1_single2$assay)) {
  res$proportion[res$assay == a] <-
    sum(pi1_single2$pi1[pi1_single2$assay == a] > 0.01) /
    sum(pi1_single2$assay == a)
}

##------------------------------------------------------------
## Wilcoxon tests: PFC vs non-PFC per assay
## (same as original loop)
##------------------------------------------------------------

wil_res <- list()
for (a in unique(pfc_single$assay)) {
  wil_res[[a]] <- wilcox.test(
    pfc_single$pi1Limma[pfc_single$assay == a & pfc_single$inPFC == TRUE],
    pfc_single$pi1Limma[pfc_single$assay == a & pfc_single$inPFC == FALSE],
    alternative = "greater"
  )
}

p_values <- data.frame(
  assay   = names(wil_res),
  p_value = sapply(wil_res, function(x) x$p.value)
)

# IMPORTANT: same default method as original (no explicit "method" argument)
p_values$fdr <- p.adjust(p_values$p_value)

# label for plotting
p_values$p_label <- ifelse(
  p_values$fdr < 0.001,
  "FDR < 0.001",
  paste("FDR =", sprintf("%.3f", p_values$fdr))
)

p_values <- as.data.frame(p_values)

##------------------------------------------------------------
## Pretty assay labels
##------------------------------------------------------------

pretty_assay <- function(x) {
  x <- gsub("Int", "Inhibitory Neurons", x)
  x <- gsub("MG",  "Microglia",           x)
  x <- gsub("Exc", "Excitatory Neurons",  x)
  x <- gsub("Oli", "Oligodendrocytes",    x)
  x <- gsub("OPC", "Oligodendrocyte\nprogenitor cells", x)
  x <- gsub("Ast", "Astrocytes",          x)
  x
}

p_values$assay2   <- pretty_assay(p_values$assay)
pfc_single$assay2 <- pretty_assay(pfc_single$assay)

##------------------------------------------------------------
## Plot
##------------------------------------------------------------

g <- ggplot(pfc_single, aes(pi1Limma, fill = inPFC)) +
  geom_density(alpha = 0.4) +
  labs(
    x    = expression(pi[1]),
    y    = "Density",
    fill = "In PFC"
  ) +
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "red"),
    breaks = c("TRUE", "FALSE")
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~assay2, scales = "free", nrow = 2) +
  geom_text(
    data = p_values,
    aes(x = Inf, y = Inf, label = p_label),
    hjust = 1.1, vjust = 2,
    size = 10,
    inherit.aes = FALSE
  ) +
  theme(
    strip.text       = element_text(face = "bold", size = 28),
    axis.title       = element_text(size = 32),
    axis.text        = element_text(size = 25),
    legend.title     = element_text(size = 27),
    legend.text      = element_text(size = 25),
    strip.background = element_rect(fill = "gray90")
  )

ggsave("pi1_pfc_single_cell_smri.pdf",g, width = 17, height = 10)
