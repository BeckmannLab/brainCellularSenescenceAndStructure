ml R/4.4.1
R

library(data.table)
library(dreamlet)
library(qvalue)
library(dplyr)
library(viridis)

myres = readRDS("single_DE_results_6.13.24.RDS")

single = myres
single2<- single[order(single$adj.P.Val, decreasing = FALSE), ]
single2 = single2[which(single2$adj.P.Val <= 0.05),]
subsetdf = as.data.frame(table(single2$feature, single2$assay))
subsetdf2 = subsetdf[order(subsetdf$Freq),]
subsetdf3=subsetdf2[which(subsetdf2$Freq >= 1),]
subsetdf3$Freq = as.numeric(subsetdf3$Freq)
subsetdf3 = subsetdf3[order(subsetdf3$Var1, decreasing = FALSE), ]
subsetdf3$Var1 = gsub("scale", "", subsetdf3$Var1)
rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
subsetdf3 <- subsetdf3[!subsetdf3$Var1 %in% rows_to_remove, ]
subsetdf3 = subsetdf3[which(subsetdf3$Var2 != "NonNeu"),]
subsetdf3 <- subsetdf3 %>%
  mutate(type = case_when(
    grepl("area", Var1, ignore.case = TRUE) ~ "area",
    grepl("thickness", Var1, ignore.case = TRUE) ~ "thickness",
    TRUE ~ "volume"
  ))

subsetdf3$Var1 = gsub("_", " ", subsetdf3$Var1)



subsetdf3$new_var = subsetdf3$Var1

subsetdf3$new_var = gsub("rh", "RH", subsetdf3$new_var)
subsetdf3$new_var = gsub("lh", "LH", subsetdf3$new_var)


subsetdf3$new_var <- gsub("LH parsorbitalis", "LH Pars Orbitalis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH parstriangularis", "RH Pars Triangularis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Putamen", "RH Putamen", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH superiorfrontal", "LH Superior Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH parsorbitalis", "RH Pars Orbitalis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Pallidum", "RH Pallidum", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH inferiorparietal", "LH Inferior Parietal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", subsetdf3$new_var)
subsetdf3$new_var <- gsub("SubCortGray", "Subcortical Gray Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH transversetemporal", "RH Transverse Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH paracentral", "RH Paracentral", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Putamen", "LH Putamen", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH entoRHinal", "LH Entorhinal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Amygdala", "LH Amygdala", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH inferiorparietal", "RH Inferior Parietal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH insula", "RH Insula", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH superiorparietal", "RH Superior Parietal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain Stem", "Brainstem", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Pallidum", "LH Pallidum", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Thalamus Proper", "RH Thalamus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH precuneus", "RH Precuneus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH middletemporal", "RH Middle Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Thalamus Proper", "LH Thalamus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH superiorfrontal", "RH Superior Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left VentralDC", "LH Ventral Diencephalon", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Amygdala", "RH Amygdala", subsetdf3$new_var)
subsetdf3$new_var <- gsub("TotalGray", "Total Gray Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH frontalpole", "RH Frontal Pole", subsetdf3$new_var)
subsetdf3$new_var <- gsub("BrainSeg", "Brain Segmentation", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH bankssts", "LH Banks STS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right vessel", "RH Vessel", subsetdf3$new_var)
subsetdf3$new_var <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH pericalcarine", "RH Pericalcarine", subsetdf3$new_var)
subsetdf3$new_var <- gsub("BrainSegNotVent", "Brain Segmentation NV", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH inferiortemporal", "LH Inferior Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right VentralDC", "RH Ventral Diencephalon", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH cuneus", "RH Cuneus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("SupraTentorial", "Supratentorial", subsetdf3$new_var)
subsetdf3$new_var <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH frontalpole", "LH Frontal Pole", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Hippocampus", "RH Hippocampus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH postcentral", "RH Postcentral", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH bankssts", "RH Banks STS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH parahippocampal", "RH Parahippocampal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH precentral", "LH Precentral", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH parahippocampal", "LH Parahippocampal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH inferiortemporal", "RH Inferior Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH temporalpole", "LH Temporal Pole", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Fifth Ventricle", "Fifth Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CSF", "Cerebrospinal Fluid", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CC Central", "Corpus Callosum Central", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH MeanThickness", "RH Mean Thickness", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH temporalpole", "RH Temporal Pole", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH supramarginal", "LH Supramarginal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH pericalcarine", "LH Pericalcarine", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CerebralWhiteMatter", "Cerebral White Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH cuneus", "LH Cuneus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH parstriangularis", "LH Pars Triangularis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH lingual", "RH Lingual", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH parsopercularis", "LH Pars Opercularis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH precuneus", "LH Precuneus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH superiortemporal", "LH Superior Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH superiorparietal", "LH Superior Parietal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH superiortemporal", "RH Superior Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH middletemporal", "LH Middle Temporal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH precentral", "RH Precentral", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH lingual", "LH Lingual", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH parsopercularis", "RH Pars Opercularis", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH supramarginal", "RH Supramarginal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH paracentral", "LH Paracentral", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Fourth Ventricle", "Fourth Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Third Ventricle", "Third Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH lateraloccipital", "LH Lateral Occipital", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LHCortex", "LH Cortex", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Cortex", "Cortex", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RHCortex", "RH Cortex", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH lateraloccipital", "RH Lateral Occipital", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Accumbens", "LH Accumbens", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH entoRHinal", "RH Entorhinal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CC Posterior", "Corpus Callosum Posterior", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right Accumbens", "RH Accumbens", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Right choroid plexus", "RH Choroid Plexus", subsetdf3$new_var)
subsetdf3$new_var <- gsub("eTIV", "Estimated Total Intracranial Volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", subsetdf3$new_var)

subsetdf3$new_var <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", subsetdf3$new_var)
subsetdf3$new_var <- gsub("SupratentorialNotVent", "Supratentorial NV", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Estimated Total Intracranial VolumeVol scaled", "Estimated Total Intracranial Volume scaled", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH CortexVol", "RH Cortex Volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("CortexVol", "Cortex Volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH CortexVol", "LH Cortex Volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", subsetdf3$new_var)

subsetdf3$new_var <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH MeanThickness thickness", "LH Mean thickness", subsetdf3$new_var)
subsetdf3$new_var <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH Mean Thickness thickness", "RH Mean thickness", subsetdf3$new_var)
subsetdf3$new_var <- gsub("RH Mean Thickness thickness", "RH Mean thickness", subsetdf3$new_var)



##sort by total
custom_breaks <- function(x) {
  max_x <- ceiling(max(x))
  c(0, seq(20, max_x, by = 20))
}

# First, determine a consistent order for cell types based on their overall frequency
cell_type_order <- subsetdf3 %>%
  group_by(Var2) %>%
  summarize(total_freq = sum(Freq)) %>%
  arrange(desc(total_freq)) %>%
  pull(Var2)

# Prepare the data
subsetdf3 <- subsetdf3 %>%
  mutate(Var2 = factor(Var2, levels = cell_type_order))

subsetdf3$celltype = subsetdf3$Var2
subsetdf3$celltype = gsub("Exc","Excitatory Neurons",subsetdf3$celltype)
subsetdf3$celltype = gsub("Int","Inhibitory Neurons",subsetdf3$celltype)
subsetdf3$celltype = gsub("MG","Microglia",subsetdf3$celltype)
subsetdf3$celltype = gsub("Oli","Oligodendrocytes",subsetdf3$celltype)
subsetdf3$celltype = gsub("OPC","Oligodendrocyte progenitor cells",subsetdf3$celltype)
subsetdf3$celltype = gsub("Ast","Astrocytes",subsetdf3$celltype)


# Create the plot
g = ggplot(subsetdf3, aes(x = Freq, y = reorder(new_var, Freq, sum))) +
  geom_segment(aes(xend = 0, yend = reorder(new_var, Freq, sum)), color = "gray") +
  geom_point(aes(color = celltype), size = 4, alpha = 0.5) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  scale_x_log10() +  # log-transformed x-axis
  theme_bw() +
  theme(
    legend.position = c(0.97, 0.02),
    legend.justification = c(1, 0),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    strip.text = element_text(size = 13),
    axis.ticks = element_line(size = 1)
  ) +
  labs(x = "Number of DEGs", y = "sMRI feature", color = "Cell Type") +
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 4)))


ggsave("degs_single_cell_smri_revised.pdf", g, width = 12, height = 8)
