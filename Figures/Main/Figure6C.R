##setup
rm(list=ls())
library(data.table)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(seriation) 
library(dplyr)
library(stringr)
library(tidyr)

##data
myres = readRDS("cor_smri_aucell_liv.RDS")
myres <- myres[assay!= "noneu"]
myres$assay = gsub("mg", "MG", myres$assay)
myres$assay = gsub("exc", "Exc", myres$assay)
myres$assay = gsub("ast", "Ast", myres$assay)
myres$assay = gsub("int", "Inh", myres$assay)
myres$assay = gsub("opc", "OPC", myres$assay)
myres$assay = gsub("oli", "Oli", myres$assay)
myres$feature <- gsub("scale", "", myres$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
myres <- myres[!myres$feature %in% rows_to_remove, ]

pm = readRDS("cor_upto5_smri_aucell.RDS")
pm <- pm[assay!="vas"]
pm$assay = gsub("micro", "MG", pm$assay)
pm$assay = gsub("exc", "Exc", pm$assay)
pm$assay = gsub("Astro", "Ast", pm$assay)
pm$assay = gsub("int", "Inh", pm$assay)
pm$assay = gsub("OPC", "OPC", pm$assay)
pm$assay = gsub("Oligo", "Oli", pm$assay)
pm$feature <- gsub("scale", "", pm$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
pm <- pm[!pm$feature %in% rows_to_remove, ]

par3 = readRDS("cor_5andup_smri_aucell.RDS")
par3 <- par3[assay!="vas"]
par3$assay = gsub("micro", "MG", par3$assay)
par3$assay = gsub("exc", "Exc", par3$assay)
par3$assay = gsub("Astro", "Ast", par3$assay)
par3$assay = gsub("int", "Inh", par3$assay)
par3$assay = gsub("OPC", "OPC", par3$assay)
par3$assay = gsub("Oligo", "Oli", par3$assay)
par3$feature <- gsub("scale", "", par3$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
par3 <- par3[!par3$feature %in% rows_to_remove, ]

myres2 = myres[order(abs(myres$spearmans), decreasing = TRUE),]
myres3 = myres2[which(myres2$pi1_sen_de > 0.01 & myres2$pi1_sc_de > 0.01),]
myres3 <- myres3 %>%
  mutate(new_column = str_extract(feature, "(?<=_)[^_]+$")) %>%
  mutate(new_column = ifelse(new_column %in% c("thickness", "area"), new_column, "volume"))

myres4 = myres3[which(myres3$p.value <= 0.05),]

myres4$new_feature = myres4$feature


myres4$new_feature = gsub("lh", "LH", myres4$new_feature)
myres4$new_feature = gsub("rh", "RH", myres4$new_feature)

myres4$new_column = gsub("volume", "Volume", myres4$new_column)
myres4$new_column = gsub("area", "Area", myres4$new_column)
myres4$new_column = gsub("thickness", "Thickness", myres4$new_column)


myres4$new_feature <- gsub("thickness", "",myres4$new_feature  )
myres4$new_feature <- gsub("area", "",myres4$new_feature  )
myres4$new_feature <- gsub("Vol", "",myres4$new_feature  )
myres4$new_feature <- gsub("_", " ",myres4$new_feature  )
myres4$new_feature <- gsub("Volume", "",myres4$new_feature  )

myres4 <- myres4 %>%
  mutate(feature = factor(feature, levels = unique(feature[order(spearmans, decreasing = TRUE)])))

# Fix naming inconsistencies
# Fix naming inconsistencies with title case replacements
myres4$new_feature <- gsub("LH parsorbitalis", "LH Pars Orbitalis", myres4$new_feature)
myres4$new_feature <- gsub("RH parstriangularis", "RH Pars Triangularis", myres4$new_feature)
myres4$new_feature <- gsub("Right Putamen", "RH Putamen", myres4$new_feature)
myres4$new_feature <- gsub("LH superiorfrontal", "LH Superior Frontal", myres4$new_feature)
myres4$new_feature <- gsub("RH parsorbitalis", "RH Pars Orbitalis", myres4$new_feature)
myres4$new_feature <- gsub("Right Pallidum", "RH Pallidum", myres4$new_feature)
myres4$new_feature <- gsub("LH inferiorparietal", "LH Inferior Parietal", myres4$new_feature)
myres4$new_feature <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", myres4$new_feature)
myres4$new_feature <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", myres4$new_feature)
myres4$new_feature <- gsub("SubCortGray", "Subcortical Gray Matter", myres4$new_feature)
myres4$new_feature <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("RH transversetemporal", "RH Transverse Temporal", myres4$new_feature)
myres4$new_feature <- gsub("RH paracentral", "RH Paracentral", myres4$new_feature)
myres4$new_feature <- gsub("Left Putamen", "LH Putamen", myres4$new_feature)
myres4$new_feature <- gsub("LH entoRHinal", "LH Entorhinal", myres4$new_feature)
myres4$new_feature <- gsub("Left Amygdala", "LH Amygdala", myres4$new_feature)
myres4$new_feature <- gsub("RH inferiorparietal", "RH Inferior Parietal", myres4$new_feature)
myres4$new_feature <- gsub("RH insula", "RH Insula", myres4$new_feature)
myres4$new_feature <- gsub("RH superiorparietal", "RH Superior Parietal", myres4$new_feature)
myres4$new_feature <- gsub("Brain Stem", "Brainstem", myres4$new_feature)
myres4$new_feature <- gsub("Left Pallidum", "LH Pallidum", myres4$new_feature)
myres4$new_feature <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", myres4$new_feature)
myres4$new_feature <- gsub("Right Thalamus Proper", "RH Thalamus", myres4$new_feature)
myres4$new_feature <- gsub("RH precuneus", "RH Precuneus", myres4$new_feature)
myres4$new_feature <- gsub("RH middletemporal", "RH Middle Temporal", myres4$new_feature)
myres4$new_feature <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", myres4$new_feature)
myres4$new_feature <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", myres4$new_feature)
myres4$new_feature <- gsub("Left Thalamus Proper", "LH Thalamus", myres4$new_feature)
myres4$new_feature <- gsub("RH superiorfrontal", "RH Superior Frontal", myres4$new_feature)
myres4$new_feature <- gsub("Left VentralDC", "LH Ventral Diencephalon", myres4$new_feature)
myres4$new_feature <- gsub("Right Amygdala", "RH Amygdala", myres4$new_feature)
myres4$new_feature <- gsub("TotalGray", "Total Gray Matter", myres4$new_feature)
myres4$new_feature <- gsub("RH frontalpole", "RH Frontal Pole", myres4$new_feature)
myres4$new_feature <- gsub("BrainSeg", "Brain Segmentation", myres4$new_feature)
myres4$new_feature <- gsub("LH bankssts", "LH Banks STS", myres4$new_feature)
myres4$new_feature <- gsub("Right vessel", "RH Vessel", myres4$new_feature)
myres4$new_feature <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", myres4$new_feature)
myres4$new_feature <- gsub("RH pericalcarine", "RH Pericalcarine", myres4$new_feature)
myres4$new_feature <- gsub("BrainSegNotVent", "Brain Segmentation NV", myres4$new_feature)
myres4$new_feature <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("LH inferiortemporal", "LH Inferior Temporal", myres4$new_feature)
myres4$new_feature <- gsub("Right VentralDC", "RH Ventral Diencephalon", myres4$new_feature)
myres4$new_feature <- gsub("RH cuneus", "RH Cuneus", myres4$new_feature)
myres4$new_feature <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("SupraTentorial", "Supratentorial", myres4$new_feature)
myres4$new_feature <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", myres4$new_feature)
myres4$new_feature <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("LH frontalpole", "LH Frontal Pole", myres4$new_feature)
myres4$new_feature <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("Right Hippocampus", "RH Hippocampus", myres4$new_feature)
myres4$new_feature <- gsub("RH postcentral", "RH Postcentral", myres4$new_feature)
myres4$new_feature <- gsub("RH bankssts", "RH Banks STS", myres4$new_feature)
myres4$new_feature <- gsub("RH parahippocampal", "RH Parahippocampal", myres4$new_feature)
myres4$new_feature <- gsub("LH precentral", "LH Precentral", myres4$new_feature)
myres4$new_feature <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", myres4$new_feature)
myres4$new_feature <- gsub("LH parahippocampal", "LH Parahippocampal", myres4$new_feature)
myres4$new_feature <- gsub("RH inferiortemporal", "RH Inferior Temporal", myres4$new_feature)
myres4$new_feature <- gsub("LH temporalpole", "LH Temporal Pole", myres4$new_feature)
myres4$new_feature <- gsub("Fifth Ventricle", "Fifth Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("CSF", "Cerebrospinal Fluid", myres4$new_feature)
myres4$new_feature <- gsub("CC Central", "Corpus Callosum Central", myres4$new_feature)
myres4$new_feature <- gsub("RH MeanThickness", "RH Mean Thickness", myres4$new_feature)
myres4$new_feature <- gsub("RH temporalpole", "RH Temporal Pole", myres4$new_feature)
myres4$new_feature <- gsub("LH supramarginal", "LH Supramarginal", myres4$new_feature)
myres4$new_feature <- gsub("LH pericalcarine", "LH Pericalcarine", myres4$new_feature)
myres4$new_feature <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", myres4$new_feature)
myres4$new_feature <- gsub("CerebralWhiteMatter", "Cerebral White Matter", myres4$new_feature)
myres4$new_feature <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("LH cuneus", "LH Cuneus", myres4$new_feature)
myres4$new_feature <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", myres4$new_feature)
myres4$new_feature <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", myres4$new_feature)
myres4$new_feature <- gsub("LH parstriangularis", "LH Pars Triangularis", myres4$new_feature)
myres4$new_feature <- gsub("RH lingual", "RH Lingual", myres4$new_feature)
myres4$new_feature <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", myres4$new_feature)
myres4$new_feature <- gsub("LH parsopercularis", "LH Pars Opercularis", myres4$new_feature)
myres4$new_feature <- gsub("LH precuneus", "LH Precuneus", myres4$new_feature)
myres4$new_feature <- gsub("LH superiortemporal", "LH Superior Temporal", myres4$new_feature)
myres4$new_feature <- gsub("LH superiorparietal", "LH Superior Parietal", myres4$new_feature)
myres4$new_feature <- gsub("RH superiortemporal", "RH Superior Temporal", myres4$new_feature)
myres4$new_feature <- gsub("LH middletemporal", "LH Middle Temporal", myres4$new_feature)
myres4$new_feature <- gsub("RH precentral", "RH Precentral", myres4$new_feature)
myres4$new_feature <- gsub("LH lingual", "LH Lingual", myres4$new_feature)
myres4$new_feature <- gsub("RH parsopercularis", "RH Pars Opercularis", myres4$new_feature)
myres4$new_feature <- gsub("RH supramarginal", "RH Supramarginal", myres4$new_feature)
myres4$new_feature <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", myres4$new_feature)
myres4$new_feature <- gsub("LH paracentral", "LH Paracentral", myres4$new_feature)
myres4$new_feature <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", myres4$new_feature)
myres4$new_feature <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", myres4$new_feature)
myres4$new_feature <- gsub("Fourth Ventricle", "Fourth Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", myres4$new_feature)
myres4$new_feature <- gsub("Third Ventricle", "Third Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("LH lateraloccipital", "LH Lateral Occipital", myres4$new_feature)
myres4$new_feature <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("LHCortex", "LH Cortex", myres4$new_feature)
myres4$new_feature <- gsub("Cortex", "Cortex", myres4$new_feature)
myres4$new_feature <- gsub("RHCortex", "RH Cortex", myres4$new_feature)
myres4$new_feature <- gsub("RH lateraloccipital", "RH Lateral Occipital", myres4$new_feature)
myres4$new_feature <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", myres4$new_feature)
myres4$new_feature <- gsub("Left Accumbens", "LH Accumbens", myres4$new_feature)
myres4$new_feature <- gsub("RH entoRHinal", "RH Entorhinal", myres4$new_feature)
myres4$new_feature <- gsub("CC Posterior", "Corpus Callosum Posterior", myres4$new_feature)
myres4$new_feature <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", myres4$new_feature)
myres4$new_feature <- gsub("Right Accumbens", "RH Accumbens", myres4$new_feature)
myres4$new_feature <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", myres4$new_feature)
myres4$new_feature <- gsub("Right choroid plexus", "RH Choroid Plexus", myres4$new_feature)
myres4$new_feature <- gsub("eTIV", "Estimated Total Intracranial Volume", myres4$new_feature)
myres4$new_feature <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", myres4$new_feature)

myres4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", myres4$new_feature)
myres4$new_feature <- gsub("SupratentorialNotVent", "Supratentorial NV", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", myres4$new_feature)
myres4$new_feature <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", myres4$new_feature)
myres4$new_feature <- gsub("Estimated Total Intracranial VolumeVol scaled", "Estimated Total Intracranial Volume scaled", myres4$new_feature)
myres4$new_feature <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", myres4$new_feature)
myres4$new_feature <- gsub("RH CortexVol", "RH Cortex Volume", myres4$new_feature)
myres4$new_feature <- gsub("CortexVol", "Cortex Volume", myres4$new_feature)
myres4$new_feature <- gsub("RH CortexVol", "LH Cortex Volume", myres4$new_feature)
myres4$new_feature <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", myres4$new_feature)
myres4$new_feature <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", myres4$new_feature)

myres4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", myres4$new_feature)
myres4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", myres4$new_feature)
myres4$new_feature <- gsub("LH MeanThickness thickness", "LH Mean thickness", myres4$new_feature)
myres4$new_feature <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", myres4$new_feature)
myres4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", myres4$new_feature)
myres4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", myres4$new_feature)

myres4$new_feature <- gsub("LHCerebral White", "LH Cerebral White Matter", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", myres4$new_feature)
myres4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", myres4$new_feature)
myres4$new_feature <- gsub("Left vessel", "LH Vessel", myres4$new_feature)


myres4$type = "Biopsy"

pm2 = pm[order(abs(pm$spearmans), decreasing = TRUE),]
pm3 = pm2[which(pm2$pi1_sen_de > 0.01 & pm2$pi1_sc_de > 0.01),]
pm3 <- pm3 %>%
  mutate(new_column = str_extract(feature, "(?<=_)[^_]+$")) %>%
  mutate(new_column = ifelse(new_column %in% c("thickness", "area"), new_column, "volume"))
pm4 = pm3[which(pm3$p.value <= 0.05),]


pm4$new_feature = pm4$feature

pm4$new_feature = gsub("lh", "LH", pm4$new_feature)
pm4$new_feature = gsub("rh", "RH", pm4$new_feature)

pm4$new_column = gsub("volume", "Volume", pm4$new_column)
pm4$new_column = gsub("area", "Area", pm4$new_column)
pm4$new_column = gsub("thickness", "Thickness", pm4$new_column)


pm4$new_feature <- gsub("thickness", "",pm4$new_feature  )
pm4$new_feature <- gsub("area", "",pm4$new_feature  )
pm4$new_feature <- gsub("Vol", "",pm4$new_feature  )
pm4$new_feature <- gsub("_", " ",pm4$new_feature  )
pm4$new_feature <- gsub("Volume", "",pm4$new_feature  )

pm4 <- pm4 %>%
  mutate(feature = factor(feature, levels = unique(feature[order(spearmans, decreasing = TRUE)])))

# Fix naming inconsistencies
pm4$new_feature <- gsub("LH parsorbitalis", "LH Pars Orbitalis", pm4$new_feature)
pm4$new_feature <- gsub("RH parstriangularis", "RH Pars Triangularis", pm4$new_feature)
pm4$new_feature <- gsub("Right Putamen", "RH Putamen", pm4$new_feature)
pm4$new_feature <- gsub("LH superiorfrontal", "LH Superior Frontal", pm4$new_feature)
pm4$new_feature <- gsub("RH parsorbitalis", "RH Pars Orbitalis", pm4$new_feature)
pm4$new_feature <- gsub("Right Pallidum", "RH Pallidum", pm4$new_feature)
pm4$new_feature <- gsub("LH inferiorparietal", "LH Inferior Parietal", pm4$new_feature)
pm4$new_feature <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", pm4$new_feature)
pm4$new_feature <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", pm4$new_feature)
pm4$new_feature <- gsub("SubCortGray", "Subcortical Gray Matter", pm4$new_feature)
pm4$new_feature <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("RH transversetemporal", "RH Transverse Temporal", pm4$new_feature)
pm4$new_feature <- gsub("RH paracentral", "RH Paracentral", pm4$new_feature)
pm4$new_feature <- gsub("Left Putamen", "LH Putamen", pm4$new_feature)
pm4$new_feature <- gsub("LH entoRHinal", "LH Entorhinal", pm4$new_feature)
pm4$new_feature <- gsub("Left Amygdala", "LH Amygdala", pm4$new_feature)
pm4$new_feature <- gsub("RH inferiorparietal", "RH Inferior Parietal", pm4$new_feature)
pm4$new_feature <- gsub("RH insula", "RH Insula", pm4$new_feature)
pm4$new_feature <- gsub("RH superiorparietal", "RH Superior Parietal", pm4$new_feature)
pm4$new_feature <- gsub("Brain Stem", "Brainstem", pm4$new_feature)
pm4$new_feature <- gsub("Left Pallidum", "LH Pallidum", pm4$new_feature)
pm4$new_feature <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", pm4$new_feature)
pm4$new_feature <- gsub("Right Thalamus Proper", "RH Thalamus", pm4$new_feature)
pm4$new_feature <- gsub("RH precuneus", "RH Precuneus", pm4$new_feature)
pm4$new_feature <- gsub("RH middletemporal", "RH Middle Temporal", pm4$new_feature)
pm4$new_feature <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", pm4$new_feature)
pm4$new_feature <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", pm4$new_feature)
pm4$new_feature <- gsub("Left Thalamus Proper", "LH Thalamus", pm4$new_feature)
pm4$new_feature <- gsub("RH superiorfrontal", "RH Superior Frontal", pm4$new_feature)
pm4$new_feature <- gsub("Left VentralDC", "LH Ventral Diencephalon", pm4$new_feature)
pm4$new_feature <- gsub("Right Amygdala", "RH Amygdala", pm4$new_feature)
pm4$new_feature <- gsub("TotalGray", "Total Gray Matter", pm4$new_feature)
pm4$new_feature <- gsub("RH frontalpole", "RH Frontal Pole", pm4$new_feature)
pm4$new_feature <- gsub("BrainSeg", "Brain Segmentation", pm4$new_feature)
pm4$new_feature <- gsub("LH bankssts", "LH Banks STS", pm4$new_feature)
pm4$new_feature <- gsub("Right vessel", "RH Vessel", pm4$new_feature)
pm4$new_feature <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", pm4$new_feature)
pm4$new_feature <- gsub("RH pericalcarine", "RH Pericalcarine", pm4$new_feature)
pm4$new_feature <- gsub("BrainSegNotVent", "Brain Segmentation NV", pm4$new_feature)
pm4$new_feature <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("LH inferiortemporal", "LH Inferior Temporal", pm4$new_feature)
pm4$new_feature <- gsub("Right VentralDC", "RH Ventral Diencephalon", pm4$new_feature)
pm4$new_feature <- gsub("RH cuneus", "RH Cuneus", pm4$new_feature)
pm4$new_feature <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("SupraTentorial", "Supratentorial", pm4$new_feature)
pm4$new_feature <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", pm4$new_feature)
pm4$new_feature <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("LH frontalpole", "LH Frontal Pole", pm4$new_feature)
pm4$new_feature <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("Right Hippocampus", "RH Hippocampus", pm4$new_feature)
pm4$new_feature <- gsub("RH postcentral", "RH Postcentral", pm4$new_feature)
pm4$new_feature <- gsub("RH bankssts", "RH Banks STS", pm4$new_feature)
pm4$new_feature <- gsub("RH parahippocampal", "RH Parahippocampal", pm4$new_feature)
pm4$new_feature <- gsub("LH precentral", "LH Precentral", pm4$new_feature)
pm4$new_feature <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", pm4$new_feature)
pm4$new_feature <- gsub("LH parahippocampal", "LH Parahippocampal", pm4$new_feature)
pm4$new_feature <- gsub("RH inferiortemporal", "RH Inferior Temporal", pm4$new_feature)
pm4$new_feature <- gsub("LH temporalpole", "LH Temporal Pole", pm4$new_feature)
pm4$new_feature <- gsub("Fifth Ventricle", "Fifth Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("CSF", "Cerebrospinal Fluid", pm4$new_feature)
pm4$new_feature <- gsub("CC Central", "Corpus Callosum Central", pm4$new_feature)
pm4$new_feature <- gsub("RH MeanThickness", "RH Mean Thickness", pm4$new_feature)
pm4$new_feature <- gsub("RH temporalpole", "RH Temporal Pole", pm4$new_feature)
pm4$new_feature <- gsub("LH supramarginal", "LH Supramarginal", pm4$new_feature)
pm4$new_feature <- gsub("LH pericalcarine", "LH Pericalcarine", pm4$new_feature)
pm4$new_feature <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", pm4$new_feature)
pm4$new_feature <- gsub("CerebralWhiteMatter", "Cerebral White Matter", pm4$new_feature)
pm4$new_feature <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("LH cuneus", "LH Cuneus", pm4$new_feature)
pm4$new_feature <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", pm4$new_feature)
pm4$new_feature <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", pm4$new_feature)
pm4$new_feature <- gsub("LH parstriangularis", "LH Pars Triangularis", pm4$new_feature)
pm4$new_feature <- gsub("RH lingual", "RH Lingual", pm4$new_feature)
pm4$new_feature <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", pm4$new_feature)
pm4$new_feature <- gsub("LH parsopercularis", "LH Pars Opercularis", pm4$new_feature)
pm4$new_feature <- gsub("LH precuneus", "LH Precuneus", pm4$new_feature)
pm4$new_feature <- gsub("LH superiortemporal", "LH Superior Temporal", pm4$new_feature)
pm4$new_feature <- gsub("LH superiorparietal", "LH Superior Parietal", pm4$new_feature)
pm4$new_feature <- gsub("RH superiortemporal", "RH Superior Temporal", pm4$new_feature)
pm4$new_feature <- gsub("LH middletemporal", "LH Middle Temporal", pm4$new_feature)
pm4$new_feature <- gsub("RH precentral", "RH Precentral", pm4$new_feature)
pm4$new_feature <- gsub("LH lingual", "LH Lingual", pm4$new_feature)
pm4$new_feature <- gsub("RH parsopercularis", "RH Pars Opercularis", pm4$new_feature)
pm4$new_feature <- gsub("RH supramarginal", "RH Supramarginal", pm4$new_feature)
pm4$new_feature <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", pm4$new_feature)
pm4$new_feature <- gsub("LH paracentral", "LH Paracentral", pm4$new_feature)
pm4$new_feature <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", pm4$new_feature)
pm4$new_feature <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", pm4$new_feature)
pm4$new_feature <- gsub("Fourth Ventricle", "Fourth Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", pm4$new_feature)
pm4$new_feature <- gsub("Third Ventricle", "Third Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("LH lateraloccipital", "LH Lateral Occipital", pm4$new_feature)
pm4$new_feature <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("LHCortex", "LH Cortex", pm4$new_feature)
pm4$new_feature <- gsub("Cortex", "Cortex", pm4$new_feature)
pm4$new_feature <- gsub("RHCortex", "RH Cortex", pm4$new_feature)
pm4$new_feature <- gsub("RH lateraloccipital", "RH Lateral Occipital", pm4$new_feature)
pm4$new_feature <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", pm4$new_feature)
pm4$new_feature <- gsub("Left Accumbens", "LH Accumbens", pm4$new_feature)
pm4$new_feature <- gsub("RH entoRHinal", "RH Entorhinal", pm4$new_feature)
pm4$new_feature <- gsub("CC Posterior", "Corpus Callosum Posterior", pm4$new_feature)
pm4$new_feature <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", pm4$new_feature)
pm4$new_feature <- gsub("Right Accumbens", "RH Accumbens", pm4$new_feature)
pm4$new_feature <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", pm4$new_feature)
pm4$new_feature <- gsub("Right choroid plexus", "RH Choroid Plexus", pm4$new_feature)
pm4$new_feature <- gsub("eTIV", "Estimated Total Intracranial Volume", pm4$new_feature)
pm4$new_feature <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", pm4$new_feature)

pm4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", pm4$new_feature)
pm4$new_feature <- gsub("SupratentorialNotVent", "Supratentorial NV", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", pm4$new_feature)
pm4$new_feature <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", pm4$new_feature)
pm4$new_feature <- gsub("Estimated Total Intracranial VolumeVol scaled", "Estimated Total Intracranial Volume scaled", pm4$new_feature)
pm4$new_feature <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", pm4$new_feature)
pm4$new_feature <- gsub("RH CortexVol", "RH Cortex Volume", pm4$new_feature)
pm4$new_feature <- gsub("CortexVol", "Cortex Volume", pm4$new_feature)
pm4$new_feature <- gsub("RH CortexVol", "LH Cortex Volume", pm4$new_feature)
pm4$new_feature <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", pm4$new_feature)
pm4$new_feature <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", pm4$new_feature)

pm4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", pm4$new_feature)
pm4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", pm4$new_feature)
pm4$new_feature <- gsub("LH MeanThickness thickness", "LH Mean thickness", pm4$new_feature)
pm4$new_feature <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", pm4$new_feature)
pm4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", pm4$new_feature)
pm4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", pm4$new_feature)

pm4$new_feature <- gsub("LHCerebral White", "LH Cerebral White Matter", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", pm4$new_feature)
pm4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", pm4$new_feature)
pm4$new_feature <- gsub("Left vessel", "LH Vessel", pm4$new_feature)

pm4$type = "Up to five"


par32 = par3[order(abs(par3$spearmans), decreasing = TRUE),]
par33 = par32[which(par32$pi1_sen_de > 0.01 & par32$pi1_sc_de > 0.01),]
par33 <- par33 %>%
  mutate(new_column = str_extract(feature, "(?<=_)[^_]+$")) %>%
  mutate(new_column = ifelse(new_column %in% c("thickness", "area"), new_column, "volume"))

par34 = par33[which(par33$p.value <= 0.05),]


par34$new_feature = par34$feature

par34$new_feature = gsub("lh", "LH", par34$new_feature)
par34$new_feature = gsub("rh", "RH", par34$new_feature)

par34$new_column = gsub("volume", "Volume", par34$new_column)
par34$new_column = gsub("area", "Area", par34$new_column)
par34$new_column = gsub("thickness", "Thickness", par34$new_column)


par34$new_feature <- gsub("thickness", "",par34$new_feature  )
par34$new_feature <- gsub("area", "",par34$new_feature  )
par34$new_feature <- gsub("Vol", "",par34$new_feature  )
par34$new_feature <- gsub("_", " ",par34$new_feature  )
par34$new_feature <- gsub("Volume", "",par34$new_feature  )

par34 <- par34 %>%
  mutate(feature = factor(feature, levels = unique(feature[order(spearmans, decreasing = TRUE)])))

# Fix naming inconsistencies
# Fix naming inconsistencies with title case replacements
par34$new_feature <- gsub("LH parsorbitalis", "LH Pars Orbitalis", par34$new_feature)
par34$new_feature <- gsub("RH parstriangularis", "RH Pars Triangularis", par34$new_feature)
par34$new_feature <- gsub("Right Putamen", "RH Putamen", par34$new_feature)
par34$new_feature <- gsub("LH superiorfrontal", "LH Superior Frontal", par34$new_feature)
par34$new_feature <- gsub("RH parsorbitalis", "RH Pars Orbitalis", par34$new_feature)
par34$new_feature <- gsub("Right Pallidum", "RH Pallidum", par34$new_feature)
par34$new_feature <- gsub("LH inferiorparietal", "LH Inferior Parietal", par34$new_feature)
par34$new_feature <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", par34$new_feature)
par34$new_feature <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", par34$new_feature)
par34$new_feature <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", par34$new_feature)
par34$new_feature <- gsub("SubCortGray", "Subcortical Gray Matter", par34$new_feature)
par34$new_feature <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", par34$new_feature)
par34$new_feature <- gsub("RH transversetemporal", "RH Transverse Temporal", par34$new_feature)
par34$new_feature <- gsub("RH paracentral", "RH Paracentral", par34$new_feature)
par34$new_feature <- gsub("Left Putamen", "LH Putamen", par34$new_feature)
par34$new_feature <- gsub("LH entoRHinal", "LH Entorhinal", par34$new_feature)
par34$new_feature <- gsub("Left Amygdala", "LH Amygdala", par34$new_feature)
par34$new_feature <- gsub("RH inferiorparietal", "RH Inferior Parietal", par34$new_feature)
par34$new_feature <- gsub("RH insula", "RH Insula", par34$new_feature)
par34$new_feature <- gsub("RH superiorparietal", "RH Superior Parietal", par34$new_feature)
par34$new_feature <- gsub("Brain Stem", "Brainstem", par34$new_feature)
par34$new_feature <- gsub("Left Pallidum", "LH Pallidum", par34$new_feature)
par34$new_feature <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", par34$new_feature)
par34$new_feature <- gsub("Right Thalamus Proper", "RH Thalamus", par34$new_feature)
par34$new_feature <- gsub("RH precuneus", "RH Precuneus", par34$new_feature)
par34$new_feature <- gsub("RH middletemporal", "RH Middle Temporal", par34$new_feature)
par34$new_feature <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", par34$new_feature)
par34$new_feature <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", par34$new_feature)
par34$new_feature <- gsub("Left Thalamus Proper", "LH Thalamus", par34$new_feature)
par34$new_feature <- gsub("RH superiorfrontal", "RH Superior Frontal", par34$new_feature)
par34$new_feature <- gsub("Left VentralDC", "LH Ventral Diencephalon", par34$new_feature)
par34$new_feature <- gsub("Right Amygdala", "RH Amygdala", par34$new_feature)
par34$new_feature <- gsub("TotalGray", "Total Gray Matter", par34$new_feature)
par34$new_feature <- gsub("RH frontalpole", "RH Frontal Pole", par34$new_feature)
par34$new_feature <- gsub("BrainSeg", "Brain Segmentation", par34$new_feature)
par34$new_feature <- gsub("LH bankssts", "LH Banks STS", par34$new_feature)
par34$new_feature <- gsub("Right vessel", "RH Vessel", par34$new_feature)
par34$new_feature <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", par34$new_feature)
par34$new_feature <- gsub("RH pericalcarine", "RH Pericalcarine", par34$new_feature)
par34$new_feature <- gsub("BrainSegNotVent", "Brain Segmentation NV", par34$new_feature)
par34$new_feature <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("LH inferiortemporal", "LH Inferior Temporal", par34$new_feature)
par34$new_feature <- gsub("Right VentralDC", "RH Ventral Diencephalon", par34$new_feature)
par34$new_feature <- gsub("RH cuneus", "RH Cuneus", par34$new_feature)
par34$new_feature <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("SupraTentorial", "Supratentorial", par34$new_feature)
par34$new_feature <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", par34$new_feature)
par34$new_feature <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("LH frontalpole", "LH Frontal Pole", par34$new_feature)
par34$new_feature <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("Right Hippocampus", "RH Hippocampus", par34$new_feature)
par34$new_feature <- gsub("RH postcentral", "RH Postcentral", par34$new_feature)
par34$new_feature <- gsub("RH bankssts", "RH Banks STS", par34$new_feature)
par34$new_feature <- gsub("RH parahippocampal", "RH Parahippocampal", par34$new_feature)
par34$new_feature <- gsub("LH precentral", "LH Precentral", par34$new_feature)
par34$new_feature <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", par34$new_feature)
par34$new_feature <- gsub("LH parahippocampal", "LH Parahippocampal", par34$new_feature)
par34$new_feature <- gsub("RH inferiortemporal", "RH Inferior Temporal", par34$new_feature)
par34$new_feature <- gsub("LH temporalpole", "LH Temporal Pole", par34$new_feature)
par34$new_feature <- gsub("Fifth Ventricle", "Fifth Ventricle", par34$new_feature)
par34$new_feature <- gsub("CSF", "Cerebrospinal Fluid", par34$new_feature)
par34$new_feature <- gsub("CC Central", "Corpus Callosum Central", par34$new_feature)
par34$new_feature <- gsub("RH MeanThickness", "RH Mean Thickness", par34$new_feature)
par34$new_feature <- gsub("RH temporalpole", "RH Temporal Pole", par34$new_feature)
par34$new_feature <- gsub("LH supramarginal", "LH Supramarginal", par34$new_feature)
par34$new_feature <- gsub("LH pericalcarine", "LH Pericalcarine", par34$new_feature)
par34$new_feature <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", par34$new_feature)
par34$new_feature <- gsub("CerebralWhiteMatter", "Cerebral White Matter", par34$new_feature)
par34$new_feature <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("LH cuneus", "LH Cuneus", par34$new_feature)
par34$new_feature <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", par34$new_feature)
par34$new_feature <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", par34$new_feature)
par34$new_feature <- gsub("LH parstriangularis", "LH Pars Triangularis", par34$new_feature)
par34$new_feature <- gsub("RH lingual", "RH Lingual", par34$new_feature)
par34$new_feature <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", par34$new_feature)
par34$new_feature <- gsub("LH parsopercularis", "LH Pars Opercularis", par34$new_feature)
par34$new_feature <- gsub("LH precuneus", "LH Precuneus", par34$new_feature)
par34$new_feature <- gsub("LH superiortemporal", "LH Superior Temporal", par34$new_feature)
par34$new_feature <- gsub("LH superiorparietal", "LH Superior Parietal", par34$new_feature)
par34$new_feature <- gsub("RH superiortemporal", "RH Superior Temporal", par34$new_feature)
par34$new_feature <- gsub("LH middletemporal", "LH Middle Temporal", par34$new_feature)
par34$new_feature <- gsub("RH precentral", "RH Precentral", par34$new_feature)
par34$new_feature <- gsub("LH lingual", "LH Lingual", par34$new_feature)
par34$new_feature <- gsub("RH parsopercularis", "RH Pars Opercularis", par34$new_feature)
par34$new_feature <- gsub("RH supramarginal", "RH Supramarginal", par34$new_feature)
par34$new_feature <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", par34$new_feature)
par34$new_feature <- gsub("LH paracentral", "LH Paracentral", par34$new_feature)
par34$new_feature <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", par34$new_feature)
par34$new_feature <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", par34$new_feature)
par34$new_feature <- gsub("Fourth Ventricle", "Fourth Ventricle", par34$new_feature)
par34$new_feature <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", par34$new_feature)
par34$new_feature <- gsub("Third Ventricle", "Third Ventricle", par34$new_feature)
par34$new_feature <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", par34$new_feature)
par34$new_feature <- gsub("LH lateraloccipital", "LH Lateral Occipital", par34$new_feature)
par34$new_feature <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", par34$new_feature)
par34$new_feature <- gsub("LHCortex", "LH Cortex", par34$new_feature)
par34$new_feature <- gsub("Cortex", "Cortex", par34$new_feature)
par34$new_feature <- gsub("RHCortex", "RH Cortex", par34$new_feature)
par34$new_feature <- gsub("RH lateraloccipital", "RH Lateral Occipital", par34$new_feature)
par34$new_feature <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", par34$new_feature)
par34$new_feature <- gsub("Left Accumbens", "LH Accumbens", par34$new_feature)
par34$new_feature <- gsub("RH entoRHinal", "RH Entorhinal", par34$new_feature)
par34$new_feature <- gsub("CC Posterior", "Corpus Callosum Posterior", par34$new_feature)
par34$new_feature <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", par34$new_feature)
par34$new_feature <- gsub("Right Accumbens", "RH Accumbens", par34$new_feature)
par34$new_feature <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", par34$new_feature)
par34$new_feature <- gsub("Right choroid plexus", "RH Choroid Plexus", par34$new_feature)
par34$new_feature <- gsub("eTIV", "Estimated Total Intracranial Volume", par34$new_feature)
par34$new_feature <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", par34$new_feature)

par34$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", par34$new_feature)
par34$new_feature <- gsub("SupratentorialNotVent", "Supratentorial NV", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", par34$new_feature)
par34$new_feature <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", par34$new_feature)
par34$new_feature <- gsub("Estimated Total Intracranial VolumeVol scaled", "Estimated Total Intracranial Volume scaled", par34$new_feature)
par34$new_feature <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", par34$new_feature)
par34$new_feature <- gsub("RH CortexVol", "RH Cortex Volume", par34$new_feature)
par34$new_feature <- gsub("CortexVol", "Cortex Volume", par34$new_feature)
par34$new_feature <- gsub("RH CortexVol", "LH Cortex Volume", par34$new_feature)
par34$new_feature <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", par34$new_feature)
par34$new_feature <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", par34$new_feature)

par34$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", par34$new_feature)
par34$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", par34$new_feature)
par34$new_feature <- gsub("LH MeanThickness thickness", "LH Mean thickness", par34$new_feature)
par34$new_feature <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", par34$new_feature)
par34$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", par34$new_feature)
par34$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", par34$new_feature)

par34$new_feature <- gsub("LHCerebral White", "LH Cerebral White Matter", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", par34$new_feature)
par34$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", par34$new_feature)
par34$new_feature <- gsub("Left vessel", "LH Vessel", par34$new_feature)

par34$type = "Five and up"


myres5 = rbind(myres4, pm4,par34)

myres5 = myres5[myres5$assay %in% c("MG","Exc"),]


mylist<- list()
for (i in c("Area", "Thickness", "Volume")){
  cur <- myres5[new_column==i]
  x1d <- dcast(cur, new_feature ~ assay, value.var="spearmans") 
  x1d.mtx <- as.matrix(x1d[,2:ncol(x1d)])
  x1d.mtx[is.na(x1d.mtx)] <- 0
  rownames(x1d.mtx) <- x1d$new_feature
  x1d.rowdist <- dist(x1d.mtx)
  x1d.coldist <- dist(t(x1d.mtx))
  x1d.roworder <- seriate(x1d.rowdist)
  x1d.colorder <- seriate(x1d.coldist)
  x1d.mtx <- x1d.mtx[unlist(x1d.roworder),unlist(x1d.colorder)]
  mylist[[i]] <- rownames(x1d.mtx)
}



mylist_liv <- list()
for (i in c("Area", "Thickness", "Volume")){
  cur <- myres4[new_column==i]
  x1d <- dcast(cur, new_feature ~ assay, value.var="spearmans") 
  x1d.mtx <- as.matrix(x1d[,2:ncol(x1d)])
  x1d.mtx[is.na(x1d.mtx)] <- 0
  rownames(x1d.mtx) <- x1d$new_feature
  x1d.rowdist <- dist(x1d.mtx)
  x1d.coldist <- dist(t(x1d.mtx))
  x1d.roworder <- seriate(x1d.rowdist)
  x1d.colorder <- seriate(x1d.coldist)
  x1d.mtx <- x1d.mtx[unlist(x1d.roworder),unlist(x1d.colorder)]
  mylist_liv[[i]] <- rownames(x1d.mtx)
}


mylist_pm <- list()
for (i in c("Area", "Thickness", "Volume")){
  cur <- pm4[new_column==i]
  x1d <- dcast(cur, new_feature ~ assay, value.var="spearmans") 
  x1d.mtx <- as.matrix(x1d[,2:ncol(x1d)])
  x1d.mtx[is.na(x1d.mtx)] <- 0
  rownames(x1d.mtx) <- x1d$new_feature
  x1d.rowdist <- dist(x1d.mtx)
  x1d.coldist <- dist(t(x1d.mtx))
  x1d.roworder <- seriate(x1d.rowdist)
  x1d.colorder <- seriate(x1d.coldist)
  x1d.mtx <- x1d.mtx[unlist(x1d.roworder),unlist(x1d.colorder)]
  mylist_pm[[i]] <- rownames(x1d.mtx)
}

mylist_par2 <- list()
for (i in c("Area", "Thickness", "Volume")){
  cur <- par34[new_column==i]
  x1d <- dcast(cur, new_feature ~ assay, value.var="spearmans") 
  x1d.mtx <- as.matrix(x1d[,2:ncol(x1d)])
  x1d.mtx[is.na(x1d.mtx)] <- 0
  rownames(x1d.mtx) <- x1d$new_feature
  x1d.rowdist <- dist(x1d.mtx)
  x1d.coldist <- dist(t(x1d.mtx))
  x1d.roworder <- seriate(x1d.rowdist)
  x1d.colorder <- seriate(x1d.coldist)
  x1d.mtx <- x1d.mtx[unlist(x1d.roworder),unlist(x1d.colorder)]
  mylist_par2[[i]] <- rownames(x1d.mtx)
}


pm4a <- pm4[new_column=="Area"]
pm4b <- pm4[new_column=="Thickness"]
pm4c <- pm4[new_column=="Volume"]
pm4a[,new_feature:=factor(new_feature, levels=mylist_pm$Area)]
pm4b[,new_feature:=factor(new_feature, levels=mylist_pm$Thickness)]
pm4c[,new_feature:=factor(new_feature, levels=mylist_pm$Volume)]

myres4a <- myres4[new_column=="Area"]
myres4b <- myres4[new_column=="Thickness"]
myres4c <- myres4[new_column=="Volume"]
myres4a[,new_feature:=factor(new_feature, levels=mylist_liv$Area)]
myres4b[,new_feature:=factor(new_feature, levels=mylist_liv$Thickness)]
myres4c[,new_feature:=factor(new_feature, levels=mylist_liv$Volume)]


myres5a <- myres5[new_column=="Area"]
myres5b <- myres5[new_column=="Thickness"]
myres5c <- myres5[new_column=="Volume"]
myres5a[,new_feature:=factor(new_feature, levels=mylist$Area)]
myres5b[,new_feature:=factor(new_feature, levels=mylist$Thickness)]
myres5c[,new_feature:=factor(new_feature, levels=mylist$Volume)]

par34a <- par34[new_column=="Area"]
par34b <- par34[new_column=="Thickness"]
par34c <- par34[new_column=="Volume"]
par34a[,new_feature:=factor(new_feature, levels=mylist_par2$Area)]
par34b[,new_feature:=factor(new_feature, levels=mylist_par2$Thickness)]
par34c[,new_feature:=factor(new_feature, levels=mylist_par2$Volume)]


myres5a$type <- factor(myres5a$type, levels = c("Up to five", "Five and up", "Biopsy"))

# For myres5b
myres5b$type <- factor(myres5b$type, levels = c("Up to five", "Five and up", "Biopsy"))

# For myres5c
myres5c$type <- factor(myres5c$type, levels = c("Up to five", "Five and up", "Biopsy"))



biopsy_res = myres5c[myres5c$type == "Up to five",]
biopsy_up = myres5c[myres5c$type == "Five and up",]

biopsy_res1 = biopsy_res[which(biopsy_res$assay == "Exc"),]
biopsy_res2 = biopsy_res[which(biopsy_res$assay == "MG"),]

biopsy_up1 = biopsy_up[which(biopsy_up$assay == "Exc"),]
biopsy_up2 = biopsy_up[which(biopsy_up$assay == "MG"),]



length(which(biopsy_res1$spearmans < 0))

# --- Plot 1: Area ---
p1 <- ggplot(myres5a, aes(x = type, y = new_feature, fill = spearmans)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    breaks = seq(-1, 1, by = 0.5),
    labels = scales::number_format(accuracy = 0.1),
    name = expression("Spearman's"~rho)
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = grid::unit(7, "cm"),
      barheight = grid::unit(0.45, "cm"),
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom"
    )
  ) +
  theme_bw() +
  ggtitle("Area") +
  theme(
    axis.text.x  = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.box.margin = margin(4, 4, 4, 4),
    legend.margin     = margin(2, 2, 2, 2),
    strip.text   = element_text(size = 14),
    plot.title   = element_text(hjust = 0.48, size = 16)
  ) +
  facet_wrap(~assay, nrow = 1, drop = FALSE) +
  labs(x = "Dataset", y = "sMRI Feature")



# --- Plot 2: Thickness ---

p2 <- ggplot(myres5b, aes(x = type, y = new_feature, fill = spearmans)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    breaks = seq(-1, 1, by = 0.5),
    labels = scales::number_format(accuracy = 0.1),
    name = expression("Spearman's"~rho)
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = grid::unit(7, "cm"),
      barheight = grid::unit(0.45, "cm"),
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom"
    )
  ) +
  theme_bw() +
  ggtitle("Thickness") +
  theme(
    axis.text.x  = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.box.margin = margin(4, 4, 4, 4),
    legend.margin     = margin(2, 2, 2, 2),
    strip.text   = element_text(size = 14),
    plot.title   = element_text(hjust = 0.48, size = 16)
  ) +
  facet_wrap(~assay, nrow = 1, drop = FALSE) +
  labs(x = "Dataset", y = "sMRI Feature")



# --- Plot 3: Volume ---
myres5c <- myres5c %>%
  group_by(new_feature) %>%
  mutate(
    mg_max = ifelse(assay == "MG", max(spearmans, na.rm = TRUE), NA)
  ) %>%
  ungroup() %>%
  mutate(
    new_feature = fct_reorder(new_feature, mg_max, .fun = max, .desc = TRUE)
  )

p3 <- ggplot(myres5c, aes(x = type, y = new_feature, fill = spearmans)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    breaks = seq(-1, 1, by = 0.5),
    labels = scales::number_format(accuracy = 0.1),
    name = expression("Spearman's"~rho)
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth  = grid::unit(8, "cm"),
      barheight = grid::unit(0.8, "cm"),
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom"
    )
  ) +
  theme_bw() +
  ggtitle("Volume") +
  theme(
    axis.text.x  = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.margin     = margin(1, 1, 1, 1),
    strip.text   = element_text(size = 14),
    plot.title   = element_text(hjust = 0.48, size = 16),
    plot.margin  = margin(5, 150, 5, 5, unit = "pt")
  ) +
  facet_wrap(~assay, nrow = 1, drop = FALSE) +
  labs(x = "Dataset", y = "sMRI Feature")



ggsave("liv_nd_uptofive_andto5_corr_smri_volume_only.pdf", p3, width = 8, height = 9)
