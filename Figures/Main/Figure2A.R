##########
# Figure 2a:
##########

ml R/4.4.1
R

library(data.table)
library(dreamlet)
library(qvalue)
library(ggplot2)

#------------------------------------------------------------
# Helper: clean / standardize sMRI feature labels
#------------------------------------------------------------
clean_feature_name <- function(x) {
  x <- gsub("rh", "RH", x)
  x <- gsub("lh", "LH", x)

  x <- gsub("LH parsorbitalis", "LH Pars Orbitalis", x)
  x <- gsub("RH parstriangularis", "RH Pars Triangularis", x)
  x <- gsub("Right Putamen", "RH Putamen", x)
  x <- gsub("LH superiorfrontal", "LH Superior Frontal", x)
  x <- gsub("RH parsorbitalis", "RH Pars Orbitalis", x)
  x <- gsub("Right Pallidum", "RH Pallidum", x)
  x <- gsub("LH inferiorparietal", "LH Inferior Parietal", x)
  x <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", x)
  x <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", x)
  x <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", x)
  x <- gsub("SubCortGray", "Subcortical Gray Matter", x)
  x <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", x)
  x <- gsub("RH transversetemporal", "RH Transverse Temporal", x)
  x <- gsub("RH paracentral", "RH Paracentral", x)
  x <- gsub("Left Putamen", "LH Putamen", x)
  x <- gsub("LH entoRHinal", "LH Entorhinal", x)
  x <- gsub("Left Amygdala", "LH Amygdala", x)
  x <- gsub("RH inferiorparietal", "RH Inferior Parietal", x)
  x <- gsub("RH insula", "RH Insula", x)
  x <- gsub("RH superiorparietal", "RH Superior Parietal", x)
  x <- gsub("Brain Stem", "Brainstem", x)
  x <- gsub("Left Pallidum", "LH Pallidum", x)
  x <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", x)
  x <- gsub("Right Thalamus Proper", "RH Thalamus", x)
  x <- gsub("RH precuneus", "RH Precuneus", x)
  x <- gsub("RH middletemporal", "RH Middle Temporal", x)
  x <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", x)
  x <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", x)
  x <- gsub("Left Thalamus Proper", "LH Thalamus", x)
  x <- gsub("RH superiorfrontal", "RH Superior Frontal", x)
  x <- gsub("Left VentralDC", "LH Ventral Diencephalon", x)
  x <- gsub("Right Amygdala", "RH Amygdala", x)
  x <- gsub("TotalGray", "Total Gray Matter", x)
  x <- gsub("RH frontalpole", "RH Frontal Pole", x)
  x <- gsub("BrainSeg", "Brain Segmentation", x)
  x <- gsub("LH bankssts", "LH Banks STS", x)
  x <- gsub("Right vessel", "RH Vessel", x)
  x <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", x)
  x <- gsub("RH pericalcarine", "RH Pericalcarine", x)
  x <- gsub("BrainSegNotVent", "Brain Segmentation NV", x)
  x <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", x)
  x <- gsub("LH inferiortemporal", "LH Inferior Temporal", x)
  x <- gsub("Right VentralDC", "RH Ventral Diencephalon", x)
  x <- gsub("RH cuneus", "RH Cuneus", x)
  x <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", x)
  x <- gsub("SupraTentorial", "Supratentorial", x)
  x <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", x)
  x <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", x)
  x <- gsub("LH frontalpole", "LH Frontal Pole", x)
  x <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", x)
  x <- gsub("Right Hippocampus", "RH Hippocampus", x)
  x <- gsub("RH postcentral", "RH Postcentral", x)
  x <- gsub("RH bankssts", "RH Banks STS", x)
  x <- gsub("RH parahippocampal", "RH Parahippocampal", x)
  x <- gsub("LH precentral", "LH Precentral", x)
  x <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", x)
  x <- gsub("LH parahippocampal", "LH Parahippocampal", x)
  x <- gsub("RH inferiortemporal", "RH Inferior Temporal", x)
  x <- gsub("LH temporalpole", "LH Temporal Pole", x)
  x <- gsub("Fifth Ventricle", "Fifth Ventricle", x)
  x <- gsub("CSF", "Cerebrospinal Fluid", x)
  x <- gsub("CC Central", "Corpus Callosum Central", x)
  x <- gsub("RH MeanThickness", "RH Mean Thickness", x)
  x <- gsub("RH temporalpole", "RH Temporal Pole", x)
  x <- gsub("LH supramarginal", "LH Supramarginal", x)
  x <- gsub("LH pericalcarine", "LH Pericalcarine", x)
  x <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", x)
  x <- gsub("CerebralWhiteMatter", "Cerebral White Matter", x)
  x <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", x)
  x <- gsub("LH cuneus", "LH Cuneus", x)
  x <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", x)
  x <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", x)
  x <- gsub("LH parstriangularis", "LH Pars Triangularis", x)
  x <- gsub("RH lingual", "RH Lingual", x)
  x <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", x)
  x <- gsub("LH parsopercularis", "LH Pars Opercularis", x)
  x <- gsub("LH precuneus", "LH Precuneus", x)
  x <- gsub("LH superiortemporal", "LH Superior Temporal", x)
  x <- gsub("LH superiorparietal", "LH Superior Parietal", x)
  x <- gsub("RH superiortemporal", "RH Superior Temporal", x)
  x <- gsub("LH middletemporal", "LH Middle Temporal", x)
  x <- gsub("RH precentral", "RH Precentral", x)
  x <- gsub("LH lingual", "LH Lingual", x)
  x <- gsub("RH parsopercularis", "RH Pars Opercularis", x)
  x <- gsub("RH supramarginal", "RH Supramarginal", x)
  x <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", x)
  x <- gsub("LH paracentral", "LH Paracentral", x)
  x <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", x)
  x <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", x)
  x <- gsub("Fourth Ventricle", "Fourth Ventricle", x)
  x <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", x)
  x <- gsub("Third Ventricle", "Third Ventricle", x)
  x <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", x)
  x <- gsub("LH lateraloccipital", "LH Lateral Occipital", x)
  x <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", x)
  x <- gsub("LHCortex", "LH Cortex", x)
  x <- gsub("Cortex", "Cortex", x)
  x <- gsub("RHCortex", "RH Cortex", x)
  x <- gsub("RH lateraloccipital", "RH Lateral Occipital", x)
  x <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", x)
  x <- gsub("Left Accumbens", "LH Accumbens", x)
  x <- gsub("RH entoRHinal", "RH Entorhinal", x)
  x <- gsub("CC Posterior", "Corpus Callosum Posterior", x)
  x <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", x)
  x <- gsub("Right Accumbens", "RH Accumbens", x)
  x <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", x)
  x <- gsub("Right choroid plexus", "RH Choroid Plexus", x)
  x <- gsub("eTIV", "Estimated Total Intracranial Volume", x)
  x <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", x)

  x <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", x)
  x <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", x)
  x <- gsub("SupratentorialNotVent", "Supratentorial NV", x)
  x <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", x)
  x <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", x)
  x <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", x)
  x <- gsub(
    "Estimated Total Intracranial VolumeVol scaled",
    "Estimated Total Intracranial Volume scaled",
    x
  )
  x <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", x)
  x <- gsub("RH CortexVol", "RH Cortex Volume", x)
  x <- gsub("CortexVol", "Cortex Volume", x)
  x <- gsub("RH CortexVol", "LH Cortex Volume", x)
  x <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", x)
  x <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", x)
  x <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", x)
  x <- gsub("LH MeanThickness thickness", "LH Mean thickness", x)
  x <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", x)
  x <- gsub("RH Mean Thickness thickness", "RH Mean thickness", x)
  x <- gsub("RH Mean Thickness thickness", "RH Mean thickness", x)

  x
}

#------------------------------------------------------------
# Read data and build feature x assay table
#------------------------------------------------------------
myres  <- readRDS("single_DE_results_6.13.24.RDS")
single <- myres

# keep significant only
single2 <- single[single$adj.P.Val <= 0.05, ]
single2 <- single2[order(single2$adj.P.Val), ]

# table of feature x assay
subsetdf <- as.data.frame(table(single2$feature, single2$assay))
colnames(subsetdf) <- c("Var1", "Var2", "Freq")

# sort and basic filtering
subsetdf2 <- subsetdf[order(subsetdf$Freq), ]
subsetdf3 <- subsetdf2[subsetdf2$Freq >= 1, ]
subsetdf3$Freq <- as.numeric(subsetdf3$Freq)
subsetdf3 <- subsetdf3[order(subsetdf3$Var1), ]

# remove "scale" tag and drop unwanted rows
subsetdf3$Var1 <- gsub("scale", "", subsetdf3$Var1)

rows_to_remove <- c(
  "Left_WM_hypointensities", "Right_WM_hypointensities",
  "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
  "WM_hypointensities", "non_WM_hypointensities",
  "MaskVol", "BrainSegVol_to_eTIV",
  "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox",
  "lhSurfaceHoles", "rhSurfaceHoles"
)

subsetdf3 <- subsetdf3[!(subsetdf3$Var1 %in% rows_to_remove), ]
subsetdf3 <- subsetdf3[subsetdf3$Var2 != "NonNeu", ]

# derive type (area / thickness / volume)
subsetdf3$type <- "volume"
idx_area      <- grepl("area", subsetdf3$Var1, ignore.case = TRUE)
idx_thickness <- grepl("thickness", subsetdf3$Var1, ignore.case = TRUE)
subsetdf3$type[idx_area]      <- "area"
subsetdf3$type[idx_thickness] <- "thickness"

# replace underscores with spaces
subsetdf3$Var1 <- gsub("_", " ", subsetdf3$Var1)

# cleaned label
subsetdf3$new_var <- clean_feature_name(subsetdf3$Var1)

#------------------------------------------------------------
# Order cell types (assays) by total frequency
#------------------------------------------------------------
cell_type_totals <- tapply(subsetdf3$Freq, subsetdf3$Var2, sum)
cell_type_order  <- names(sort(cell_type_totals, decreasing = TRUE))

subsetdf3$Var2 <- factor(subsetdf3$Var2, levels = cell_type_order)

#------------------------------------------------------------
# Nicely labelled cell type variable
#------------------------------------------------------------
subsetdf3$celltype <- as.character(subsetdf3$Var2)
subsetdf3$celltype <- gsub("Exc", "Excitatory Neurons", subsetdf3$celltype)
subsetdf3$celltype <- gsub("Int", "Inhibitory Neurons", subsetdf3$celltype)
subsetdf3$celltype <- gsub("MG",  "Microglia", subsetdf3$celltype)
subsetdf3$celltype <- gsub("Oli", "Oligodendrocytes", subsetdf3$celltype)
subsetdf3$celltype <- gsub(
  "OPC", "Oligodendrocyte progenitor cells", subsetdf3$celltype
)
subsetdf3$celltype <- gsub("Ast", "Astrocytes", subsetdf3$celltype)

# keep features significant in at least 2 cell types
tab_feat_ct <- table(subsetdf3$new_var, subsetdf3$celltype)
qualified   <- names(which(rowSums(tab_feat_ct) >= 2))
filtered_df <- subsetdf3[subsetdf3$new_var %in% qualified, ]

#------------------------------------------------------------
# Plot
#------------------------------------------------------------
g <- ggplot(filtered_df,
            aes(x = Freq,
                y = reorder(new_var, Freq, sum))) +
  geom_segment(aes(xend = 0,
                   yend = reorder(new_var, Freq, sum)),
               color = "gray") +
  geom_point(aes(color = celltype), size = 6, alpha = 0.5) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  scale_x_log10() +
  theme_bw() +
  theme(
    legend.position      = c(0.97, 0.02),
    legend.justification = c(1, 0),
    legend.box.just      = "right",
    legend.margin        = margin(6, 6, 6, 6),
    legend.background    = element_rect(fill = "white", color = "black"),
    legend.key.size      = unit(1.5, "lines"),
    legend.text          = element_text(size = 25),
    legend.title         = element_text(size = 27),
    axis.text            = element_text(size = 25),
    axis.title           = element_text(size = 28),
    strip.text           = element_text(size = 20),
    axis.ticks           = element_line(size = 1)
  ) +
  labs(
    x     = "Number of DEGs",
    y     = "sMRI feature",
    color = "Cell Type"
  ) +
  guides(
    alpha = "none",
    color = guide_legend(override.aes = list(size = 6))
  )

ggsave(
  "degs_single_cell_smri_revised.pdf",
  g, width = 13, height = 6
)
