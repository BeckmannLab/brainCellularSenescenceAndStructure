library(qvalue)
library(dplyr)
library(ggplot2)
library(limma)
library(data.table)

resFinal=fread("allSTRDEGresults_DE_noODC_bulk_01.10.25.txt",data.table=FALSE) 
myres = resFinal
rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")
rows_to_remove2 = intersect(rows_to_remove, myres$feature)

# Remove rows with specified row names
df2 <- myres[!myres$feature %in% rows_to_remove2, ]

summary = df2 %>% 
	group_by(feature) %>%
	summarise(
		nonSignif=sum(adj.P.Val>0.05),
		signif=sum(adj.P.Val<=0.05),
		pi1=1-propTrueNull(P.Value)
		) %>%
	arrange(desc(signif))

summary %>% arrange(desc(pi1))

summary=as.data.frame(summary)

length(which(summary$pi1 > 0))
pfc <- c( "caudalmiddlefrontal", "lateralorbitofrontal", "medialorbitofrontal", 
           "parsopercularis", "parsorbitalis", "parstriangularis", "rostralmiddlefrontal",
           "superiorfrontal", "frontalpole")
library(rex)
pfc_regex=rex(or(pfc))

summary[grep(pfc_regex,summary$feature,perl=TRUE),] %>% arrange(desc(pi1))
summary$inPFC=FALSE
summary$inPFC[grep(pfc_regex,summary$feature,perl=TRUE)]=TRUE

wil_res = wilcox.test(summary$pi1[summary$inPFC == TRUE],summary$pi1[summary$inPFC == FALSE], alternative = "greater")

p_values <- wil_res$p.value

g = ggplot(summary, aes(pi1, fill = inPFC)) +
  geom_density(alpha = 0.4) +
  labs(x = expression(pi[1]), y = "Density", fill = "In PFC") +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                    breaks = c("TRUE", "FALSE")) +
  annotate("text", x = Inf, y = Inf, label = "p-value = 0.27",
           hjust = 1.1, vjust = 3.5, size = 7) +  # slightly lower
  theme_bw() + 
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 14),   
    axis.text.x = element_text(size = 14),   
    axis.title.y = element_text(size = 16),  
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
ggsave("pi1_pfc_bulk_smri.pdf",g, width = 8 , height = 8)

