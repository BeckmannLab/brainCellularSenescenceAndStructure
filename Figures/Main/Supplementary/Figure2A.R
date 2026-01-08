# Load needed packages
library(data.table)
library(ggplot2)
library(grid)  

## -----------------------------
## 1. Load and prepare RNA-seq counts
## -----------------------------

genedata_org <- fread(
  "geneCounts_2022-07-26.txt",
  data.table = FALSE
)

genedata <- as.data.frame(genedata_org[, 2:ncol(genedata_org)])
rownames(genedata) <- genedata_org[, 1]

# Keep samples with sufficient total counts
genedata <- genedata[, colSums(genedata) > 1.7e7]

## -----------------------------
## 2. Load and prepare imaging + covariates
## -----------------------------

str_data_summarized <- readRDS(
  "str_data_summarized_JAN172024.RDS"
)
str_data_summarized <- str_data_summarized$E

cov_str_summarized <- readRDS(
  "cov_str_summarized_JAN172024.RDS"
)

# Align covariates to imaging data
cov_str_summarized <- cov_str_summarized[
  match(colnames(str_data_summarized), cov_str_summarized$subject_id),
]

# Sample info
info_all <- fread(
  "info_all_withImaging_2022-07-26.txt",
  sep = "\t", header = TRUE, data.table = FALSE
)

pheTable <- fread(
  "lbpPheFile.txt",
  data.table = FALSE
)

info_all <- merge(
  info_all,
  pheTable,
  by.x = "IID_ISMMS",
  by.y = "iid",
  suffixes = c("", "All")
)

# Order info to match RNA samples
info_all2 <- info_all[match(colnames(genedata), info_all$SAMPLE_ISMMS), ]

# Merge imaging covariates and sample info
fullInfo <- merge(
  cov_str_summarized,
  info_all2,
  by.x   = "subject_id",
  by.y   = "IID_ISMMS",
  all.y  = TRUE,
  suffixes = c(".imaging", ".RNA")
)

fullInfo2 <- fullInfo[match(colnames(genedata), fullInfo$SAMPLE_ISMMS), ]

# Keep subjects with both RNA and imaging
fullInfo2 <- fullInfo2[fullInfo2$subject_id %in% colnames(str_data_summarized), ]

# Align RNA and imaging to the same set/order of subjects
genedata <- genedata[, fullInfo2$SAMPLE_ISMMS]
str_data_summarized <- str_data_summarized[, fullInfo2$subject_id]

# Add imaging features as columns
fullInfo2 <- cbind(fullInfo2, t(str_data_summarized))

## -----------------------------
## 3. Build data for chosen feature & gene
## -----------------------------

# Feature and gene you want to plot
feature_name <- "rh_WhiteSurfArea_area"
gene_id      <- "ENSG00000224272.2"

subset <- fullInfo2[, c(feature_name, "SAMPLE_ISMMS")]
rownames(subset) <- subset$SAMPLE_ISMMS
colnames(subset)[1] <- "feature"

# Take that one gene across samples
genedata_subset <- as.data.frame(
  t(genedata[rownames(genedata) == gene_id, ])
)

# Merge gene expression with imaging feature
all <- merge(genedata_subset, subset, by = "row.names")
colnames(all)[2] <- "gene"  # rename gene column to something simple

## -----------------------------
## 4. Plot and save
## -----------------------------

g <- ggplot(all, aes(x = gene, y = feature)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(
    x = "AC114730.3",
    y = "Right Hemisphere White Surface Area"
  ) +
  theme(
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 18),
    axis.ticks.length = unit(0.1, "cm")
  )

ggsave(
  "example_DE.pdf",
  g, width = 9, height = 8
)
