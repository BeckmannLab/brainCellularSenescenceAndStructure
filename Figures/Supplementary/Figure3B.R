#################
library(dplyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

myres = readRDS("new_annotation_bulk_pi1_sig_GO_by_features.RDS")

myres$feature <- sub("_[^_]*$", "", myres$feature)

rows_to_remove <- c("Left_WM_hypointensities", "Right_WM_hypointensities",
                     "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
                     "WM_hypointensities", "non_WM_hypointensities",
                     "MaskVol", "BrainSegVol_to_eTIV",
                     "MaskVol_to_eTIV", "SurfaceHoles", "SupraTentorialVolNotVentVox", "lhSurfaceHoles", "rhSurfaceHoles")



# Remove rows with specified row names
myres <- myres[!myres$feature %in% rows_to_remove, ]


myres$significant = ifelse(myres$BH <= 0.05, TRUE, FALSE)
myres2 = myres[which(myres$BH <= 0.05),]
myres2 = as.data.frame(myres2)
myres2$feature <- sub("_[^_]*$", "", myres2$feature)
myres2 = myres2[myres2$celltype != "NonNeu",]
sig_GO_infl= myres2[grep("inflammatory response", myres2$Term),]
sig_GO_stress= myres2[grep("stress", myres2$Term),]
sig_GO_oxidative_phosphorylation= myres2[grep("oxidative phosphorylation", myres2$Term),]

sig_GO_misc = myres2[grep("GO:0090398|GO:0001302|GO:0001300|GO:0044838|GO:0007050|GO:0034640|GO:0034641|GO:0006281|GO:0006974|GO:0000723|GO:0032200|GO:0006979|GO:0016049|GO:0008285|GO:0051301|GO:0007569|GO:0042127|GO:0008283|GO:0051726|GO:0007568",myres2$GO.ID),]


single_sen = rbind(sig_GO_infl, sig_GO_stress, sig_GO_oxidative_phosphorylation, sig_GO_misc)
myres2 = single_sen

res = data.frame(GO.ID = unique(myres2$GO.ID), mg= NA, exc= NA, oli = NA, OPC = NA, int = NA, ast = NA)
for(x in unique(myres2$GO.ID)){
  sub = myres2[myres2$GO.ID == x,]
  tab= table(sub$celltype)
  res$mg[res$GO.ID == x] = tab["MG"]
  res$exc[res$GO.ID == x] = tab["Exc"]
  res$oli[res$GO.ID == x] = tab["Oli"]
  res$OPC[res$GO.ID == x] = tab["OPC"]
  res$int[res$GO.ID == x] = tab["Int"]
  res$ast[res$GO.ID == x] = tab["Ast"]
}

res$ncell_type = rowSums(res[,2:7]/res[,2:7], na.rm = TRUE)
res$mean_go_term = rowMeans(res[,2:7],na.rm = TRUE)


res2 = merge(unique(myres2[,c("GO.ID", "Term")]), res)

res2 = res2[order(res2$mean_go_term),]
res2 = res2[order(res2$ncell_type),]

res2[is.na(res2)] <- 0

df_long <- res2 %>%
  pivot_longer(cols = c(mg, exc, oli, OPC, int, ast), 
               names_to = "CellType", 
               values_to = "Count")

# Calculate total count for each term
df_long <- df_long %>%
  group_by(Term) %>%
  mutate(Total = sum(Count)) %>%
  ungroup()

df_long[df_long == 0] <- NA

df_long$CellType = gsub("mg", "Microglia", df_long$CellType)
df_long$CellType = gsub("exc", "Excitatory Neurons", df_long$CellType)
df_long$CellType = gsub("oli", "Oligodendrocytes", df_long$CellType)
df_long$CellType = gsub("int", "Inhibitory Neurons", df_long$CellType)
df_long$CellType = gsub("ast", "Astrocytes", df_long$CellType)
df_long$CellType = gsub("OPC", "Oligodendrocyte progenitor cells", df_long$CellType)

color_scheme <- c(
  "Microglia" = "#1b9e77",
  "Excitatory Neurons" = "#e7298a",
  "Oligodendrocytes" = "#66a61e",
  "Inhibitory Neurons" = "#d95f02",
  "Astrocytes" = "#7570b3",
  "Oligodendrocyte progenitor cells" = "#e6ab02"
)

df_long[which(df_long$Term == "stress-activated protein kinase signalin..."),]
df_long$Term = gsub("regulation of cell population proliferat...","regulation of cell population proliferation",df_long$Term)
df_long$Term = gsub("negative regulation of cell population p...","negative regulation of cell population proliferation",df_long$Term)
df_long$Term = gsub("stress-activated protein kinase signalin...","stress-activated protein kinase signaling cascade",df_long$Term)

# Create the stacked bar plot
b = ggplot(df_long, aes(x = Count, y = reorder(Term, Total), fill = CellType)) +
  geom_bar(stat = "identity") +
 geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5),
            color = "black", size = 6) +
  labs(x = "Total Count", y = "GO Term", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  scale_fill_manual(values = color_scheme) +
  # scale_x_continuous(expand = c(0, 0), 
  #                    breaks = seq(0, max(df_long$Total), by = 20)) +  # Added more ticks
  coord_cartesian(clip = "off")
  
ggsave("single_GO_sen.pdf", b, height = 10, width = 15,limitsize = FALSE)
