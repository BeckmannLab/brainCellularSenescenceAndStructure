
# setup 
  library(dreamlet)
  library(data.table)
  library(qvalue)
  library(ggplot2)
  library(ggthemes)
  library(forcats)
  library(dplyr)
  library(scales)

# pi1
pi1_single = readRDS("single_pi1_per_feature_celltype_09.16.24.RDS")

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")

pi1_single <- pi1_single[!pi1_single$feature %in% rows_to_remove, ]
pi1_single = pi1_single[which(pi1_single$assay != "NonNeu"),]
pi1_single2 = pi1_single

# pi1_single2 = pi1_single[which(pi1_single$pi1Limma > 0.05),]
# pi1_single2 <- pi1_single2[pi1_single2$pi1Limma != 0, ]
pi1_single2 <- pi1_single2 %>%
  mutate(type = case_when(
    grepl("area", feature, ignore.case = TRUE) ~ "area",
    grepl("thickness", feature, ignore.case = TRUE) ~ "thickness",
    TRUE ~ "volume"
  ))

pfc <- data.frame(feature = c("caudalmiddlefrontal", "lateralorbitofrontal", "medialorbitofrontal","parsopercularis", "parsorbitalis", "parstriangularis", "rostralmiddlefrontal","superiorfrontal", "frontalpole"),feature1 = NA, feature2 = NA)

for (feature in pfc$feature){
	tmp = grep(feature, pi1_single2$feature, value = T)
	pfc$feature1[pfc$feature == feature] = tmp[1]
	pfc$feature2[pfc$feature == feature] = tmp[2]
}

pfc2 = data.frame(feature = c(pfc$feature1, pfc$feature2), inPFC = TRUE)
pfc_single = merge(pfc2, pi1_single2, by = "feature", all.y = TRUE)
pfc_single$inPFC[which(is.na(pfc_single$inPFC))] = FALSE

res = data.frame(assay= unique(pi1_single2$assay), proportion = NA)
for (i in unique(pi1_single2$assay)){
  res$proportion[res$assay==i]=sum(pi1_single2$pi1[pi1_single2$assay==i]>0.01)/sum(pi1_single2$assay==i)
}

wil_res = list()
for (i in unique(pfc_single$assay)){
  wil_res[[i]] = wilcox.test(pfc_single$pi1Limma[pfc_single$assay == i & pfc_single$inPFC == TRUE],pfc_single$pi1Limma[pfc_single$assay == i & pfc_single$inPFC == FALSE], alternative = "greater")
}

p_values <- data.frame(
  assay = names(wil_res),
  p_value = sapply(wil_res, function(x) x$p.value)
)

p_values$fdr = p.adjust(p_values$p_value)
p_values <- p_values %>%
  mutate(p_label = ifelse(fdr < 0.001, 
                          "FDR < 0.001",
                          paste("FDR =", sprintf("%.3f", fdr))))

p_values = as.data.frame(p_values)

p_values$assay2 = p_values$assay
p_values$assay2 = gsub("Int", "Inhibitory Neurons",p_values$assay2)
p_values$assay2 = gsub("MG", "Microglia",p_values$assay2)
p_values$assay2 = gsub("Exc", "Excitatory Neurons",p_values$assay2)
p_values$assay2 = gsub("Oli", "Oligodendrocytes",p_values$assay2)
p_values$assay2 = gsub("OPC", "Oligodendrocyte\nprogenitor cells",p_values$assay2)
p_values$assay2 = gsub("Ast", "Astrocytes",p_values$assay2)

pfc_single$assay2= pfc_single$assay
pfc_single$assay2 = gsub("Int", "Inhibitory Neurons",pfc_single$assay2)
pfc_single$assay2 = gsub("MG", "Microglia",pfc_single$assay2)
pfc_single$assay2 = gsub("Exc", "Excitatory Neurons",pfc_single$assay2)
pfc_single$assay2 = gsub("Oli", "Oligodendrocytes",pfc_single$assay2)
pfc_single$assay2 = gsub("OPC", "Oligodendrocyte\nprogenitor cells",pfc_single$assay2)
pfc_single$assay2 = gsub("Ast", "Astrocytes",pfc_single$assay2)


g <- ggplot(pfc_single, aes(pi1Limma, fill = inPFC)) +
  geom_density(alpha = 0.4) +
  labs(
    x = expression(pi[1]),
    y = "Density",
    fill = "In PFC"
  ) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +  # 2 decimal points on x-axis
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "red"),
    breaks = c("TRUE", "FALSE")
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~assay2, scales = "fixed", nrow = 3) +
  geom_text(
    data = p_values, 
    aes(x = Inf, y = Inf, label = p_label),
    hjust = 1.1, vjust = 2,
    size = 10,
    inherit.aes = FALSE
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 28),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 27),
    legend.text = element_text(size = 25),
    strip.background = element_rect(fill = "gray90")
  ) 
  
ggsave("pi1_pfc_single_cell_smri_same_axes.pdf",g, width = 13 , height = 12)
