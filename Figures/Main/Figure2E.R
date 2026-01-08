##############################
# Fig 2E
##############################

library(ggplot2)

##---------------------------
## Load data
##---------------------------

bulk <- readRDS("bulk_with_no_sex_sn_summ_9.20.24.RDS")
prot <- readRDS("prot_sn_summ_9.20.24.RDS")

bulk <- as.data.frame(bulk)
prot <- as.data.frame(prot)

##---------------------------
## Filter & prepare bulk
##---------------------------

bulk_sub <- bulk[
  bulk$pi1.sc   >= 0.01 &
  bulk$pi1.bulk >= 0.01 &
  bulk$fdr      <= 0.05,
]

# order by |rho| desc, then feature
bulk_sub <- bulk_sub[order(abs(bulk_sub$spearmans), decreasing = TRUE), ]
bulk_sub <- bulk_sub[order(bulk_sub$feature, decreasing = TRUE), ]

bulk_sub$omic <- "Bulk RNAseq"
bulk_sub <- bulk_sub[, c("feature", "assay", "spearmans", "omic")]

##---------------------------
## Filter & prepare proteomics
##---------------------------

prot_sub <- prot[
  prot$pi1.sc  >= 0.01 &
  prot$pi1.prt >= 0.01 &
  prot$fdr     <= 0.05,
]

prot_sub <- prot_sub[order(abs(prot_sub$spearmans), decreasing = TRUE), ]
prot_sub <- prot_sub[order(prot_sub$feature, decreasing = TRUE), ]

prot_sub$omic <- "MS"
prot_sub <- prot_sub[, c("feature", "assay", "spearmans", "omic")]

##---------------------------
## Combine & clean features
##---------------------------

all <- rbind(bulk_sub, prot_sub)

rows_to_remove <- c(
  "Left_WM_hypointensities", "Right_WM_hypointensities",
  "Left_non_WM_hypointensities", "Right_non_WM_hypointensities",
  "WM_hypointensities", "non_WM_hypointensities",
  "MaskVol", "BrainSegVol_to_eTIV",
  "MaskVol_to_eTIV", "SupraTentorialVolNotVentVox",
  "lhSurfaceHoles", "rhSurfaceHoles", "SurfaceHoles"
)

all <- all[!(all$feature %in% rows_to_remove), ]
all <- all[all$assay != "NonNeu", ]

##---------------------------
## Derived variables
##---------------------------

all$abs_spearman <- abs(all$spearmans)

all$celltype <- all$assay
all$assay    <- gsub("Int", "Inh", all$assay)

# Optional color scheme (not used in current plot, but kept for reference)
color_scheme <- c(
  "MG"  = "#1b9e77",
  "Exc" = "#e7298a",
  "Oli" = "#66a61e",
  "Inh" = "#d95f02",
  "Ast" = "#7570b3",
  "OPC" = "#e6ab02"
)

##---------------------------
## Plot
##---------------------------

all$omic <- gsub("MS", "Mass spectrometry", all$omic)

b <- ggplot(all, aes(x = assay, y = abs_spearman)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(
    x = "Cell Type",
    y = expression("Absolute correlation with snRNAseq (Spearman's " * rho * ")")
  ) +
  theme_bw() +
  facet_wrap(~omic) +
  theme(
    axis.title  = element_text(size = 24),
    axis.text   = element_text(size = 20),
    strip.text  = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 17)
  )

ggsave("bulk_prot_sc_cor.pdf",b, width = 10, height = 10)
