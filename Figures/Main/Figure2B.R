##########
# 2b:
##########

# setup
library(dreamlet)
library(data.table)
library(qvalue)
library(ggplot2)

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
pi1_single <- readRDS("single_pi1_per_feature_celltype_09.16.24.RDS")

##------------------------------------------------------------
## Drop unwanted features / cell types
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
## Region type (area / thickness / volume)
##------------------------------------------------------------
pi1_single2$type <- "volume"
idx_area      <- grepl("area",      pi1_single2$feature, ignore.case = TRUE)
idx_thickness <- grepl("thickness", pi1_single2$feature, ignore.case = TRUE)

pi1_single2$type[idx_area]      <- "area"
pi1_single2$type[idx_thickness] <- "thickness"

##------------------------------------------------------------
## Clean feature labels
##------------------------------------------------------------
clean_feature_labels <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("rh", "Right Hemisphere", x)
  x <- gsub("lh", "Left Hemisphere", x)

  x <- gsub(
    "Right Hemisphere rostralanteriorcingulate thickness",
    "Right Hemisphere rostral anterior cingulate", x
  )
  x <- gsub("parsorbitalis",           "pars orbitalis",          x)
  x <- gsub("parstriangularis",        "pars triangularis",       x)
  x <- gsub("CerebralWhiteMatter",     "Cerebral White Matter",   x)
  x <- gsub("caudalanteriorcingulate", "caudal anterior cingulate", x)
  x <- gsub(
    "HemisphereCerebralWhiteMatter",
    "Hemisphere Cerebral White Matter", x
  )
  x <- gsub("HemisphereCerebral", "Hemisphere Cerebral", x)
  x <- gsub("MeanThickness",      "Mean Thickness",      x)

  x
}

pi1_single2$feature <- clean_feature_labels(pi1_single2$feature)

##------------------------------------------------------------
## Proportion of signatures with pi1 > 0.01 per assay
##------------------------------------------------------------
threshold <- 0.01

# tapply returns mean of logical (TRUE/FALSE) ⇒ proportion
prop_by_assay <- tapply(pi1_single2$pi1 > threshold,
                        pi1_single2$assay,
                        mean)

res <- data.frame(
  assay     = names(prop_by_assay),
  proportion = as.numeric(prop_by_assay),
  row.names = NULL
)

# rename inhibitory label
res$assay <- gsub("Int", "Inh", res$assay)

##------------------------------------------------------------
## Bar plot
##------------------------------------------------------------
g <- ggplot(res, aes(x = reorder(assay, -proportion), y = proportion)) +
  geom_bar(stat = "identity", fill = "grey40") +
  theme_bw() +
  theme(
    axis.text.x      = element_text(size = 25),
    axis.text.y      = element_text(size = 25),
    axis.title.x     = element_text(size = 28),
    axis.title.y     = element_text(size = 28),
    legend.position  = "none",
    axis.ticks       = element_line(size = 1),
    axis.ticks.length = unit(0.10, "cm")
  ) +
  labs(
    x = "Cell Type",
    y = expression(paste("Proportion of signatures with ", pi[1], " > 0.01"))
  )

ggsave("pi1_single_cell_smri.pdf",g, width = 10, height = 8)

