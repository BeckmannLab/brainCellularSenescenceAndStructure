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

pm = readRDS("corr_pmAucellscore_smri_features.RDS")
pm <- pm[assay!="nonneu"]
pm$assay = gsub("mg", "MG", pm$assay)
pm$assay = gsub("exc", "Exc", pm$assay)
pm$assay = gsub("ast", "Ast", pm$assay)
pm$assay = gsub("int", "Inh", pm$assay)
pm$assay = gsub("opc", "OPC", pm$assay)
pm$assay = gsub("oli", "Oli", pm$assay)

pm$feature <- gsub("scale", "", pm$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
pm <- pm[!pm$feature %in% rows_to_remove, ]


bs = readRDS("corr_psychencodeAucellscore_smri_features.RDS")
bs$assay = bs$celltype
bs <- bs[assay!= "noneu"]
bs$assay = gsub("Micro", "MG", bs$assay)
# bs$assay = gsub("exc", "Exc", bs$assay)
bs$assay = gsub("Astro", "Ast", bs$assay)
bs$assay = gsub("Int", "Inh", bs$assay)
# bs$assay = gsub("opc", "OPC", bs$assay)
bs$assay = gsub("Oligo", "Oli", bs$assay)

bs$feature <- gsub("scale", "", bs$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

# Remove rows with specified row names
bs <- bs[!bs$feature %in% rows_to_remove, ]



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

pm4$type = "Bank"


###brainscope data


bs2 = bs[order(abs(bs$spearmans), decreasing = TRUE),]
bs3 = bs2[which(bs2$pi1_sen_de > 0.01 & bs2$pi1_sc_de > 0.01),]
bs3 <- bs3 %>%
  mutate(new_column = str_extract(feature, "(?<=_)[^_]+$")) %>%
  mutate(new_column = ifelse(new_column %in% c("thickness", "area"), new_column, "volume"))
bs4 = bs3[which(bs3$p.value <= 0.05),]

bs4$new_feature = bs4$feature

bs4$new_feature = gsub("lh", "LH", bs4$new_feature)
bs4$new_feature = gsub("rh", "RH", bs4$new_feature)

bs4$new_column = gsub("volume", "Volume", bs4$new_column)
bs4$new_column = gsub("area", "Area", bs4$new_column)
bs4$new_column = gsub("thickness", "Thickness", bs4$new_column)


bs4$new_feature <- gsub("thickness", "",bs4$new_feature  )
bs4$new_feature <- gsub("area", "",bs4$new_feature  )
bs4$new_feature <- gsub("Vol", "",bs4$new_feature  )
bs4$new_feature <- gsub("_", " ",bs4$new_feature  )
bs4$new_feature <- gsub("Volume", "",bs4$new_feature  )

bs4 <- bs4 %>%
  mutate(feature = factor(feature, levels = unique(feature[order(spearmans, decreasing = TRUE)])))

# Fix naming inconsistencies
# Fix naming inconsistencies with title case replacements
bs4$new_feature <- gsub("LH parsorbitalis", "LH Pars Orbitalis", bs4$new_feature)
bs4$new_feature <- gsub("RH parstriangularis", "RH Pars Triangularis", bs4$new_feature)
bs4$new_feature <- gsub("Right Putamen", "RH Putamen", bs4$new_feature)
bs4$new_feature <- gsub("LH superiorfrontal", "LH Superior Frontal", bs4$new_feature)
bs4$new_feature <- gsub("RH parsorbitalis", "RH Pars Orbitalis", bs4$new_feature)
bs4$new_feature <- gsub("Right Pallidum", "RH Pallidum", bs4$new_feature)
bs4$new_feature <- gsub("LH inferiorparietal", "LH Inferior Parietal", bs4$new_feature)
bs4$new_feature <- gsub("RH isthmuscingulate", "RH Isthmus Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("LH caudalmiddlefrontal", "LH Caudal Middle Frontal", bs4$new_feature)
bs4$new_feature <- gsub("Right Cerebellum Cortex", "RH Cerebellum Cortex", bs4$new_feature)
bs4$new_feature <- gsub("SubCortGray", "Subcortical Gray Matter", bs4$new_feature)
bs4$new_feature <- gsub("LH isthmuscingulate", "LH Isthmus Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("RH transversetemporal", "RH Transverse Temporal", bs4$new_feature)
bs4$new_feature <- gsub("RH paracentral", "RH Paracentral", bs4$new_feature)
bs4$new_feature <- gsub("Left Putamen", "LH Putamen", bs4$new_feature)
bs4$new_feature <- gsub("LH entoRHinal", "LH Entorhinal", bs4$new_feature)
bs4$new_feature <- gsub("Left Amygdala", "LH Amygdala", bs4$new_feature)
bs4$new_feature <- gsub("RH inferiorparietal", "RH Inferior Parietal", bs4$new_feature)
bs4$new_feature <- gsub("RH insula", "RH Insula", bs4$new_feature)
bs4$new_feature <- gsub("RH superiorparietal", "RH Superior Parietal", bs4$new_feature)
bs4$new_feature <- gsub("Brain Stem", "Brainstem", bs4$new_feature)
bs4$new_feature <- gsub("Left Pallidum", "LH Pallidum", bs4$new_feature)
bs4$new_feature <- gsub("RH lateralorbitofrontal", "RH Lateral Orbitofrontal", bs4$new_feature)
bs4$new_feature <- gsub("Right Thalamus Proper", "RH Thalamus", bs4$new_feature)
bs4$new_feature <- gsub("RH precuneus", "RH Precuneus", bs4$new_feature)
bs4$new_feature <- gsub("RH middletemporal", "RH Middle Temporal", bs4$new_feature)
bs4$new_feature <- gsub("Left Cerebellum Cortex", "LH Cerebellum Cortex", bs4$new_feature)
bs4$new_feature <- gsub("LH lateralorbitofrontal", "LH Lateral Orbitofrontal", bs4$new_feature)
bs4$new_feature <- gsub("Left Thalamus Proper", "LH Thalamus", bs4$new_feature)
bs4$new_feature <- gsub("RH superiorfrontal", "RH Superior Frontal", bs4$new_feature)
bs4$new_feature <- gsub("Left VentralDC", "LH Ventral Diencephalon", bs4$new_feature)
bs4$new_feature <- gsub("Right Amygdala", "RH Amygdala", bs4$new_feature)
bs4$new_feature <- gsub("TotalGray", "Total Gray Matter", bs4$new_feature)
bs4$new_feature <- gsub("RH frontalpole", "RH Frontal Pole", bs4$new_feature)
bs4$new_feature <- gsub("BrainSeg", "Brain Segmentation", bs4$new_feature)
bs4$new_feature <- gsub("LH bankssts", "LH Banks STS", bs4$new_feature)
bs4$new_feature <- gsub("Right vessel", "RH Vessel", bs4$new_feature)
bs4$new_feature <- gsub("BrainSegNotVentSurf", "Brain Segmentation NVS", bs4$new_feature)
bs4$new_feature <- gsub("RH pericalcarine", "RH Pericalcarine", bs4$new_feature)
bs4$new_feature <- gsub("BrainSegNotVent", "Brain Segmentation NV", bs4$new_feature)
bs4$new_feature <- gsub("RH caudalanteriorcingulate", "RH Caudal Anterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("LH inferiortemporal", "LH Inferior Temporal", bs4$new_feature)
bs4$new_feature <- gsub("Right VentralDC", "RH Ventral Diencephalon", bs4$new_feature)
bs4$new_feature <- gsub("RH cuneus", "RH Cuneus", bs4$new_feature)
bs4$new_feature <- gsub("RH rostralanteriorcingulate", "RH Rostral Anterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("SupraTentorial", "Supratentorial", bs4$new_feature)
bs4$new_feature <- gsub("SupraTentorialNotVent", "Supratentorial No Ventricles", bs4$new_feature)
bs4$new_feature <- gsub("RH posteriorcingulate", "RH Posterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("LH frontalpole", "LH Frontal Pole", bs4$new_feature)
bs4$new_feature <- gsub("LH caudalanteriorcingulate", "LH Caudal Anterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("Right Hippocampus", "RH Hippocampus", bs4$new_feature)
bs4$new_feature <- gsub("RH postcentral", "RH Postcentral", bs4$new_feature)
bs4$new_feature <- gsub("RH bankssts", "RH Banks STS", bs4$new_feature)
bs4$new_feature <- gsub("RH parahippocampal", "RH Parahippocampal", bs4$new_feature)
bs4$new_feature <- gsub("LH precentral", "LH Precentral", bs4$new_feature)
bs4$new_feature <- gsub("Right Cerebellum White Matter", "RH Cerebellum White Matter", bs4$new_feature)
bs4$new_feature <- gsub("LH parahippocampal", "LH Parahippocampal", bs4$new_feature)
bs4$new_feature <- gsub("RH inferiortemporal", "RH Inferior Temporal", bs4$new_feature)
bs4$new_feature <- gsub("LH temporalpole", "LH Temporal Pole", bs4$new_feature)
bs4$new_feature <- gsub("Fifth Ventricle", "Fifth Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("CSF", "Cerebrospinal Fluid", bs4$new_feature)
bs4$new_feature <- gsub("CC Central", "Corpus Callosum Central", bs4$new_feature)
bs4$new_feature <- gsub("RH MeanThickness", "RH Mean Thickness", bs4$new_feature)
bs4$new_feature <- gsub("RH temporalpole", "RH Temporal Pole", bs4$new_feature)
bs4$new_feature <- gsub("LH supramarginal", "LH Supramarginal", bs4$new_feature)
bs4$new_feature <- gsub("LH pericalcarine", "LH Pericalcarine", bs4$new_feature)
bs4$new_feature <- gsub("RHCerebralWhiteMatter", "RH Cerebral White Matter", bs4$new_feature)
bs4$new_feature <- gsub("CerebralWhiteMatter", "Cerebral White Matter", bs4$new_feature)
bs4$new_feature <- gsub("LH rostralanteriorcingulate", "LH Rostral Anterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("LH cuneus", "LH Cuneus", bs4$new_feature)
bs4$new_feature <- gsub("RH rostralmiddlefrontal", "RH Rostral Middle Frontal", bs4$new_feature)
bs4$new_feature <- gsub("CC Mid Anterior", "Corpus Callosum Mid Anterior", bs4$new_feature)
bs4$new_feature <- gsub("LH parstriangularis", "LH Pars Triangularis", bs4$new_feature)
bs4$new_feature <- gsub("RH lingual", "RH Lingual", bs4$new_feature)
bs4$new_feature <- gsub("LHCerebralWhiteMatter", "LH Cerebral White Matter", bs4$new_feature)
bs4$new_feature <- gsub("LH parsopercularis", "LH Pars Opercularis", bs4$new_feature)
bs4$new_feature <- gsub("LH precuneus", "LH Precuneus", bs4$new_feature)
bs4$new_feature <- gsub("LH superiortemporal", "LH Superior Temporal", bs4$new_feature)
bs4$new_feature <- gsub("LH superiorparietal", "LH Superior Parietal", bs4$new_feature)
bs4$new_feature <- gsub("RH superiortemporal", "RH Superior Temporal", bs4$new_feature)
bs4$new_feature <- gsub("LH middletemporal", "LH Middle Temporal", bs4$new_feature)
bs4$new_feature <- gsub("RH precentral", "RH Precentral", bs4$new_feature)
bs4$new_feature <- gsub("LH lingual", "LH Lingual", bs4$new_feature)
bs4$new_feature <- gsub("RH parsopercularis", "RH Pars Opercularis", bs4$new_feature)
bs4$new_feature <- gsub("RH supramarginal", "RH Supramarginal", bs4$new_feature)
bs4$new_feature <- gsub("Left Cerebellum White Matter", "LH Cerebellum White Matter", bs4$new_feature)
bs4$new_feature <- gsub("LH paracentral", "LH Paracentral", bs4$new_feature)
bs4$new_feature <- gsub("LH rostralmiddlefrontal", "LH Rostral Middle Frontal", bs4$new_feature)
bs4$new_feature <- gsub("LH posteriorcingulate", "LH Posterior Cingulate", bs4$new_feature)
bs4$new_feature <- gsub("Fourth Ventricle", "Fourth Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("RH medialorbitofrontal", "RH Medial Orbitofrontal", bs4$new_feature)
bs4$new_feature <- gsub("Third Ventricle", "Third Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("Left Lateral Ventricle", "LH Lateral Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("LH lateraloccipital", "LH Lateral Occipital", bs4$new_feature)
bs4$new_feature <- gsub("Right Lateral Ventricle", "RH Lateral Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("LHCortex", "LH Cortex", bs4$new_feature)
bs4$new_feature <- gsub("Cortex", "Cortex", bs4$new_feature)
bs4$new_feature <- gsub("RHCortex", "RH Cortex", bs4$new_feature)
bs4$new_feature <- gsub("RH lateraloccipital", "RH Lateral Occipital", bs4$new_feature)
bs4$new_feature <- gsub("CC Mid Posterior", "Corpus Callosum Mid Posterior", bs4$new_feature)
bs4$new_feature <- gsub("Left Accumbens", "LH Accumbens", bs4$new_feature)
bs4$new_feature <- gsub("RH entoRHinal", "RH Entorhinal", bs4$new_feature)
bs4$new_feature <- gsub("CC Posterior", "Corpus Callosum Posterior", bs4$new_feature)
bs4$new_feature <- gsub("Left Inf Lat Vent", "LH Inferior Lateral Ventricle", bs4$new_feature)
bs4$new_feature <- gsub("Right Accumbens", "RH Accumbens", bs4$new_feature)
bs4$new_feature <- gsub("RH caudalmiddlefrontal", "RH Caudal Middle Frontal", bs4$new_feature)
bs4$new_feature <- gsub("Right choroid plexus", "RH Choroid Plexus", bs4$new_feature)
bs4$new_feature <- gsub("eTIV", "Estimated Total Intracranial Volume", bs4$new_feature)
bs4$new_feature <- gsub("EstimatedTotalIntraCranial", "Estimated Total Intracranial Volume", bs4$new_feature)

bs4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", bs4$new_feature)
bs4$new_feature <- gsub("SupratentorialNotVent", "Supratentorial NV", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationVolNotVent", "Brain Segmentation Volume NV", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationVolNotVentSurf", "Brain Segmentation Volume NVS", bs4$new_feature)
bs4$new_feature <- gsub("Total Gray MatterVol", "Total Gray Matter Volume", bs4$new_feature)
bs4$new_feature <- gsub("Estimated Total Intracranial VolumeVol scaled", "Estimated Total Intracranial Volume scaled", bs4$new_feature)
bs4$new_feature <- gsub("Brain Segmentation Volume NVSurf", "Brain Segmentation Volume NVS", bs4$new_feature)
bs4$new_feature <- gsub("RH CortexVol", "RH Cortex Volume", bs4$new_feature)
bs4$new_feature <- gsub("CortexVol", "Cortex Volume", bs4$new_feature)
bs4$new_feature <- gsub("RH CortexVol", "LH Cortex Volume", bs4$new_feature)
bs4$new_feature <- gsub("LHCerebral White MatterVol", "LH Cerebral White Matter volume", bs4$new_feature)
bs4$new_feature <- gsub("RH Cerebral White MatterVol", "RH Cerebral White Matter volume", bs4$new_feature)

bs4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", bs4$new_feature)
bs4$new_feature <- gsub("Cerebral White MatterVol", "Cerebral White Matter volume", bs4$new_feature)
bs4$new_feature <- gsub("LH MeanThickness thickness", "LH Mean thickness", bs4$new_feature)
bs4$new_feature <- gsub("LH Cerebral White Matter Volume", "LH Cerebral White Matter volume", bs4$new_feature)
bs4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", bs4$new_feature)
bs4$new_feature <- gsub("RH Mean Thickness thickness", "RH Mean thickness", bs4$new_feature)

bs4$new_feature <- gsub("LHCerebral White", "LH Cerebral White Matter", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationNotVent", "Brain Segmentation NV", bs4$new_feature)
bs4$new_feature <- gsub("Brain SegmentationNotVentSurf", "Brain Segmentation NVS", bs4$new_feature)
bs4$new_feature <- gsub("Left vessel", "LH Vessel", bs4$new_feature)

bs4$type = "brainSCOPE"

bs4$celltype = NULL
bs4$percent = NULL

myres5 = rbind(pm4, myres4,bs4)


# ---- 1) Function to get row/column orders from seriation ----
get_orders <- function(dt) {
  x <- dcast(dt, new_feature ~ assay, value.var = "spearmans")
  m <- as.matrix(x[, -1])
  m[is.na(m)] <- 0
  rownames(m) <- x$new_feature
  row_ord <- rownames(m)[unlist(seriate(dist(m)))]
  col_ord <- colnames(m)[unlist(seriate(dist(t(m))))]
  list(rows = row_ord, cols = col_ord)
}

# ---- 2) Get orders per measurement type ----
orders <- lapply(c("Area", "Thickness", "Volume"),
                 \(k) get_orders(myres5[new_column == k]))
names(orders) <- c("Area", "Thickness", "Volume")

# ---- 3) Subset data ----
myres5a <- myres5[new_column=="Area"]
myres5b <- myres5[new_column=="Thickness"]
myres5c <- myres5[new_column=="Volume"]

# ---- 4) Apply row orders ----
myres5a[, new_feature := factor(new_feature, levels = orders$Area$rows)]
myres5b[, new_feature := factor(new_feature, levels = orders$Thickness$rows)]
myres5c[, new_feature := factor(new_feature, levels = orders$Volume$rows)]

# ---- 5) Apply facet (assay) orders to match screenshot ----
thick_order  <- c("MG","Inh","Oli","Ast","OPC","Exc")     # A. Thickness
area_order   <- c("MG","Oli","Inh","Ast","OPC","Exc")     # B. Area
vol_order    <- c("OPC","MG","Inh","Ast","Exc","Oli")     # C. Volume

myres5a[, assay := factor(assay, levels = area_order)]
myres5b[, assay := factor(assay, levels = thick_order)]
myres5c[, assay := factor(assay, levels = vol_order)]

# ---- 6) Dataset order (x-axis) ----
ds_order <- c("Biopsy","Bank","brainSCOPE")
myres5a[, type := factor(type, levels = ds_order)]
myres5b[, type := factor(type, levels = ds_order)]
myres5c[, type := factor(type, levels = ds_order)]

# ---- 7) Plot function ----
plot_panel <- function(dt, title) {
  ggplot(dt, aes(x = type, y = new_feature, fill = spearmans)) +
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
        barwidth = grid::unit(7, "cm"),
        barheight = grid::unit(0.45, "cm"),
        title.position = "left",
        title.hjust = 0.5,
        label.position = "bottom"
      )
    ) +
    theme_bw() +
    ggtitle(title) +
    theme(
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 11),
      legend.box.margin = margin(4, 4, 4, 4),
      legend.margin = margin(2, 2, 2, 2),
      strip.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.48, size = 16)
    ) +
    facet_wrap(~assay, nrow = 1, drop = FALSE) +
    labs(x = "Dataset", y = "sMRI Feature")
}

# ---- 8) Generate plots ----
p1 <- plot_panel(myres5b, "Thickness") # A
p2 <- plot_panel(myres5a, "Area")      # B
p3 <- plot_panel(myres5c, "Volume")    # C

# ---- 9) Combine panels in order ----
all <- p1 | p2 | p3

all <- (p1 | p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("liv_pm_bs_corr_smri_only_volume.pdf", all, width = 8, height = 10)
