library(data.table)
library(SingleCellExperiment)
library(zellkonverter)
library(dreamlet)
library(zenith)
library(DelayedArray)
library(GSEABase)
library(cowplot)
library(tidyverse)
library(kableExtra)
library(qvalue)
library(scattermore)
library(RColorBrewer)
library(corrplot)
library(viridis)
library(dplyr)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(ggrepel)
library(Seurat)
library(assertthat)
library(Matrix)
library(foreach)
library(AUCell)
library(matrixStats)
library(tidyr)
library(broom)
library(scales)
library(magrittr)
library(patchwork)

################################################################################################################
# row 1 – AUCell boxplots
################################################################################################################

##########
# GSE115301
##########
meta <- fread("GSE115301_Growing_Sen_10x_metadata.txt")
AUCX2 <- readRDS("AUCell_GSE115301_1.16.24.RDS")
AUC_mat <- t(getAUC(AUCX2))
AUC_mat <- as.data.frame(AUC_mat)

all <- merge(AUC_mat, meta, by.x = "row.names", by.y = "V1")

all$status <- "non-senescent"
all$status[all$V2 == "Senescence1"] <- "senescent"
all$status <- factor(all$status, levels = c("senescent", "non-senescent"))
all$ds <- "GSE115301"
all_2 <- all
table(all_2$status)

# senescent non-senescent 
#       273           207 

wilcox_results <- all_2 %>%
  group_by(ds) %>%
  summarise(
    p_value = wilcox.test(sen_list ~ status, alternative = "greater", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_text = paste("p =", scientific(p_value, digits = 2)))

b <- ggplot(all_2, aes(x = status, y = sen_list, fill = ds)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  labs(
    x = "Cell State",
    y = "AUCell Score",
    title = "Teo YV et al., 2019 Cell Rep"
  ) +
  theme_bw() +
  facet_wrap(~ds) +
  scale_fill_manual(values = c("white")) +
  geom_text(
    data = wilcox_results,
    aes(
      label = p_text,
      x = max(as.numeric(all_2$status)) + 0.5,
      y = max(all_2$sen_list) * 0.999
    ),
    hjust = 1,
    vjust = 1,
    size = 10,
    inherit.aes = FALSE
  ) +
  theme(
    axis.title   = element_text(size = 24),
    axis.text    = element_text(size = 22),
    strip.text   = element_blank(),
    plot.title   = element_text(size = 26, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

##########
# GSE175533
##########
AUCX2 <- readRDS("AUCell_GSE175533_1.16.24.RDS")
meta <- readRDS("meta_Data.RDS")

AUC_mat <- t(getAUC(AUCX2))
AUC_mat <- as.data.frame(AUC_mat)

all <- merge(AUC_mat, meta, by = "row.names")

all$status <- "non-senescent"
all$status[all$exp == "WT"] <- "senescent"
all$status <- factor(all$status, levels = c("senescent", "non-senescent"))
all$ds <- "GSE175533"
all_3 <- as.data.frame(all)

table(all_3$status)

# senescent non-senescent 
#      8953          3066 

wilcox_results <- all_3 %>%
  group_by(ds) %>%
  summarise(
    p_value = wilcox.test(sen_list ~ status, alternative = "greater", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_text = paste("p =", scientific(p_value, digits = 2)))

d <- ggplot(all_3, aes(x = status, y = sen_list, fill = ds)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  labs(
    x = "Cell State",
    y = "AUCell Score",
    title = "Chan et al., 2022 Elife"
  ) +
  theme_bw() +
  facet_wrap(~ds) +
  scale_fill_manual(values = c("white")) +
  geom_text(
    data = wilcox_results,
    aes(
      label = p_text,
      x = max(as.numeric(all_3$status)) + 0.5,
      y = max(all_3$sen_list) * 0.999
    ),
    hjust = 1,
    vjust = 1,
    size = 10,
    inherit.aes = FALSE
  ) +
  theme(
    axis.title   = element_text(size = 24),
    axis.text    = element_text(size = 22),
    strip.text   = element_blank(),
    plot.title   = element_text(size = 26, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

##########
# GSE119807
##########
AUCX2 <- readRDS("AUCell_GSE119807_1.23.24.RDS")
meta <- readRDS("meta_for_counts.RDS")

AUC_mat <- t(getAUC(AUCX2))
AUC_mat <- as.data.frame(AUC_mat)

all <- merge(AUC_mat, meta, by.x = "row.names", by.y = "cell")

all$status <- "non-senescent"
all$status[all$type == "sen"] <- "senescent"
all$status <- factor(all$status, levels = c("senescent", "non-senescent"))
all$ds <- "GSE119807"
all_1 <- all

table(all_1$status)

# senescent non-senescent 
#       800           400 

wilcox_results <- all_1 %>%
  group_by(ds) %>%
  summarise(
    p_value = wilcox.test(sen_list ~ status, alternative = "greater", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_text = paste("p =", scientific(p_value, digits = 2)))

g <- ggplot(all_1, aes(x = status, y = sen_list, fill = ds)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  labs(
    x = "Cell State",
    y = "AUCell Score",
    title = "Tang et al., 2019 Protein Cell"
  ) +
  theme_bw() +
  facet_wrap(~ds) +
  scale_fill_manual(values = c("white")) +
  geom_text(
    data = wilcox_results,
    aes(
      label = p_text,
      x = max(as.numeric(all_1$status)) + 0.5,
      y = max(all_1$sen_list) * 0.99
    ),
    hjust = 1,
    vjust = 1,
    size = 10,
    inherit.aes = FALSE
  ) +
  theme(
    axis.title   = element_text(size = 24),
    axis.text    = element_text(size = 22),
    strip.text   = element_blank(),
    plot.title   = element_text(size = 26, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

n <- (b + d + g) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30))

################################################################################################################
# row 2 – odds ratios vs AUCell thresholds
################################################################################################################

##########
# GSE115301
##########
meta <- fread("GSE115301_Growing_Sen_10x_metadata.txt")
AUCX2 <- readRDS("AUCell_GSE115301_1.16.24.RDS")
V1 <- c("0.7", "0.8", "0.9", "0.95", "0.99")
V2 <- c("30percent", "20percent", "10percent", "5percent", "1percent")
thresholds <- as.data.frame(cbind(V1, V2))
thresholds$V1 <- as.numeric(thresholds$V1)

df <- data.frame(test = character(), odds_ratio = numeric(), fdr = numeric())
for (x in 1:nrow(thresholds)) {
  AUC_mat <- t(getAUC(AUCX2))
  AUC_mat <- as.data.frame(AUC_mat)
  AUC_mat$V1 <- rownames(AUC_mat)
  threshold <- quantile(AUC_mat$sen_list, thresholds[x, 1])
  AUC_mat$sen <- AUC_mat$sen_list >= threshold
  AUC_mat$percent <- thresholds[x, 1]
  AUC_mat$cell <- rownames(AUC_mat)

  all <- merge(AUC_mat, meta, by.x = "row.names", by.y = "V1")
  all$V4 <- 0
  all$V4[all$V2 == "Senescence1"] <- 1
  all$category <- "sen"
  all$category[all$V2 != "Senescence1"] <- "control"

  te <- fisher.test(table(all$sen, all$category), alternative = "greater")

  df_row <- data.frame(
    test       = thresholds[x, 2],
    odds_ratio = te$estimate,
    fdr        = te$p.value
  )
  df <- rbind(df_row, df)
}

df$test <- factor(df$test, levels = c("1percent", "5percent", "10percent", "20percent", "30percent"))
df_1 <- df
df_1$type <- "True Senescence"

predicted_sen <- readRDS("enrichment_or_12.16.25.RDS")
predicted_sen2 <- predicted_sen[, c("var2", "odds_ratio", "p_value")]
colnames(predicted_sen2) <- c("test", "odds_ratio", "fdr")
predicted_sen2$type <- "Predicted Senescence"

# clean for plotting
df_1$test <- gsub("percent", "", df_1$test)
df_1$test <- factor(df_1$test, levels = c("1", "5", "10", "20", "30"))

predicted_sen2$test <- gsub("percent", "", predicted_sen2$test)
predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))

true_sen_1 <- ggplot(df_1, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

# ensure full set of levels 1,5,10,20,30 for predicted
predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))
desired_levels <- c("1", "5", "10", "20", "30")
missing_levels <- setdiff(desired_levels, unique(predicted_sen2$test))

missing_rows <- data.frame(
  test       = missing_levels,
  odds_ratio = NA,
  fdr        = NA,
  type       = "Predicted Senescence"
)
predicted_sen2_full <- rbind(predicted_sen2, missing_rows)
predicted_sen2_full$test <- factor(predicted_sen2_full$test, levels = desired_levels)

pred_sen_1 <- ggplot(predicted_sen2_full, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio",
    title = "Teo YV et al., 2019 Cell Rep"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

##########
# GSE175533
##########
AUCX2 <- readRDS("AUCell_GSE175533_1.16.24.RDS")
meta <- readRDS("meta_Data.RDS")
V1 <- c("0.7", "0.8", "0.9", "0.95", "0.99")
V2 <- c("30percent", "20percent", "10percent", "5percent", "1percent")
thresholds <- as.data.frame(cbind(V1, V2))
thresholds$V1 <- as.numeric(thresholds$V1)

df <- data.frame(test = character(), odds_ratio = numeric(), fdr = numeric())
for (x in 1:nrow(thresholds)) {
  AUC_mat <- t(getAUC(AUCX2))
  AUC_mat <- as.data.frame(AUC_mat)
  AUC_mat$V1 <- rownames(AUC_mat)
  threshold <- quantile(AUC_mat$sen_list, thresholds[x, 1])
  AUC_mat$sen <- AUC_mat$sen_list >= threshold
  AUC_mat$percent <- thresholds[x, 1]
  AUC_mat$cell <- rownames(AUC_mat)

  all <- merge(AUC_mat, meta, by = "row.names")
  all$V4 <- 1
  all$V4[all$exp == "tert"] <- 0
  all$category <- "control"
  all$category[all$exp == "WT"] <- "senescent"

  te <- fisher.test(table(all$sen, all$category), alternative = "greater")

  df_row <- data.frame(
    test       = thresholds[x, 2],
    odds_ratio = te$estimate,
    fdr        = te$p.value
  )
  df <- rbind(df_row, df)
}

df_3 <- df
df_3$type <- "True Senescence"

predicted_sen <- readRDS("enrichment_or_12.16.25.RDS")
predicted_sen2 <- predicted_sen[, c("var2", "odds_ratio", "p_value")]
colnames(predicted_sen2) <- c("test", "odds_ratio", "fdr")
predicted_sen2$type <- "Predicted Senescence"

predicted_sen2$test <- gsub("percent", "", predicted_sen2$test)
predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))

df_3$test <- gsub("percent", "", df_3$test)
df_3$test <- factor(df_3$test, levels = c("1", "5", "10", "20", "30"))

true_sen_2 <- ggplot(df_3, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))
desired_levels <- c("1", "5", "10", "20", "30")
missing_levels <- setdiff(desired_levels, unique(predicted_sen2$test))

missing_rows <- data.frame(
  test       = missing_levels,
  odds_ratio = NA,
  fdr        = NA,
  type       = "Predicted Senescence"
)
predicted_sen2_full <- rbind(predicted_sen2, missing_rows)
predicted_sen2_full$test <- factor(predicted_sen2_full$test, levels = desired_levels)

pred_sen_2 <- ggplot(predicted_sen2, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio",
    title = "Chan et al., 2022 Elife"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

##########
# GSE119807
##########
AUCX2 <- readRDS("AUCell_GSE119807_1.23.24.RDS")
meta <- readRDS("meta_for_counts.RDS")
V1 <- c("0.7", "0.8", "0.9", "0.95", "0.99")
V2 <- c("30percent", "20percent", "10percent", "5percent", "1percent")
thresholds <- as.data.frame(cbind(V1, V2))
thresholds$V1 <- as.numeric(thresholds$V1)

df <- data.frame(test = character(), odds_ratio = numeric(), fdr = numeric())
for (x in 1:nrow(thresholds)) {
  AUC_mat <- t(getAUC(AUCX2))
  AUC_mat <- as.data.frame(AUC_mat)
  AUC_mat$V1 <- rownames(AUC_mat)
  threshold <- quantile(AUC_mat$sen_list, thresholds[x, 1])
  AUC_mat$sen <- AUC_mat$sen_list >= threshold
  AUC_mat$percent <- thresholds[x, 1]
  AUC_mat$cell <- rownames(AUC_mat)

  all <- merge(AUC_mat, meta, by.x = "row.names", by.y = "cell")
  all$V4 <- 0
  all$V4[all$type == "sen"] <- 1
  all$category <- "sen"
  all$category[all$type != "sen"] <- "control"

  te <- fisher.test(table(all$sen, all$category), alternative = "greater")

  df_row <- data.frame(
    test       = thresholds[x, 2],
    odds_ratio = te$estimate,
    fdr        = te$p.value
  )
  df <- rbind(df_row, df)
}

df$test <- factor(df$test, levels = c("1percent", "5percent", "10percent", "20percent", "30percent"))
df_2 <- df
df_2$type <- "True Senescence"

predicted_sen <- readRDS("enrichment_or_12.16.25.RDS")
predicted_sen2 <- predicted_sen[, c("var2", "odds_ratio", "p_value")]
colnames(predicted_sen2) <- c("test", "odds_ratio", "fdr")
predicted_sen2$type <- "Predicted Senescence"

df_2$test <- gsub("percent", "", df_2$test)
df_2$test <- factor(df_2$test, levels = c("1", "5", "10", "20", "30"))

predicted_sen2$test <- gsub("percent", "", predicted_sen2$test)
predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))

true_sen_3 <- ggplot(df_2, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

predicted_sen2$test <- factor(predicted_sen2$test, levels = c("1", "5", "10", "20", "30"))
desired_levels <- c("1", "5", "10", "20", "30")
missing_levels <- setdiff(desired_levels, unique(predicted_sen2$test))

missing_rows <- data.frame(
  test       = missing_levels,
  odds_ratio = NA,
  fdr        = NA,
  type       = "Predicted Senescence"
)
predicted_sen2_full <- rbind(predicted_sen2, missing_rows)
predicted_sen2_full$test <- factor(predicted_sen2_full$test, levels = desired_levels)

pred_sen_3 <- ggplot(predicted_sen2_full, aes(x = test, y = odds_ratio)) +
  geom_line(aes(group = 1), color = "black", size = 1) +
  geom_point(aes(fill = fdr > 0.05), color = "black", size = 3, shape = 21) +
  scale_fill_manual(
    values = c("black", "white"),
    labels = list(expression("FDR" <= 0.05), expression("FDR" > 0.05)),
    name   = "Significance"
  ) +
  labs(
    x = "Senescence cell proportion threshold",
    y = "Odds Ratio",
    title = "Tang et al., 2019 Protein Cell"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text      = element_text(size = 20),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 20),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.box.margin     = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  expand_limits(y = 1)

predicted_row1 <- ((pred_sen_1 + pred_sen_2 + pred_sen_3) +
  plot_layout(axes = "collect")) +
  theme(legend.position = "bottom")

true_row1 <- ((true_sen_1 + true_sen_2 + true_sen_3) +
  plot_layout(axes = "collect")) +
  theme(legend.position = "bottom")

################################################################################################################
# row 3 – DE logFC correlations (true vs predicted)
################################################################################################################

##########
# GSE115301
##########
df <- readRDS(paste0("correct_edger", "_", "true_sen", "_fit.RDS"))
df$symbol <- rownames(df)
df$sig <- df$FDR <= 0.05
true_sen <- df

var <- read.delim("var_de_aucell_mg_exc.txt", header = FALSE)
var <- var[, 1:3]
var <- unique(var)

i <- 5
file_path <- paste0("correct_edger", "_", var[i, 2], "_fit.RDS")
df <- readRDS(file_path)
df$test <- var[i, 2]
df$symbol <- rownames(df)
all2 <- df
all2$sig <- all2$FDR <= 0.05
predicted_sen <- all2

all <- merge(true_sen, predicted_sen, by = c("symbol"), all.y = TRUE)

correlations <- all %>%
  summarise(
    test_result = list(cor.test(logFC.x, logFC.y, method = "spearman", use = "complete.obs"))
  ) %>%
  mutate(
    rho    = sapply(test_result, function(x) x$estimate),
    p_value = sapply(test_result, function(x) x$p.value),
    label  = sprintf("rho = %.2f, p = %.3f", rho, p_value)
  )

b <- ggplot(all, aes(x = logFC.x, y = logFC.y)) +
  geom_point(color = "grey40", alpha = 0.7) +
  theme_bw(base_size = 14) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "True Senescence LogFC",
    y = "Predicted Senescence\nLogFC"
  ) +
  theme(
    text       = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20)
  ) +
  geom_text(
    data = correlations,
    aes(label = label),
    x = 0.95 * max(all$logFC.x, na.rm = TRUE),
    y = 0.95 * max(all$logFC.y, na.rm = TRUE),
    hjust = 1, vjust = 1,
    size = 10
  )

##########
# GSE175533
##########
var <- read.delim("var_de_aucell_mg_exc.txt", header = FALSE)
var <- var[, 1:3]
var <- unique(var)

df <- readRDS(paste0("correct_edger", "_", "true_sen", "_fit.RDS"))
df$symbol <- rownames(df)
all2 <- df
all2$sig <- all2$FDR <= 0.05
true_sen <- all2

i <- 1
file_path <- paste0("correct_edger", "_", var[i, 2], "_fit.RDS")
df <- readRDS(file_path)
df$test <- var[i, 2]
df$symbol <- rownames(df)
all2 <- df
all2$sig <- all2$FDR <= 0.05
predicted_sen <- all2

all <- merge(true_sen, predicted_sen, by = c("symbol"), all.y = TRUE)

correlations <- all %>%
  summarise(
    test_result = list(cor.test(logFC.x, logFC.y, method = "spearman", use = "complete.obs"))
  ) %>%
  mutate(
    rho    = sapply(test_result, function(x) x$estimate),
    p_value = sapply(test_result, function(x) x$p.value),
    label  = sprintf("rho = %.2f, p = %.3f", rho, p_value)
  )

h <- ggplot(all, aes(x = logFC.x, y = logFC.y)) +
  geom_point(color = "grey40", alpha = 0.7) +
  theme_bw(base_size = 14) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "True Senescence LogFC",
    y = "Predicted Senescence\nLogFC"
  ) +
  theme(
    text       = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20)
  ) +
  geom_text(
    data = correlations,
    aes(label = label),
    x = 0.95 * max(all$logFC.x, na.rm = TRUE),
    y = 0.95 * max(all$logFC.y, na.rm = TRUE),
    hjust = 1, vjust = 1,
    size = 10
  )

##########
# GSE119807
##########
df <- readRDS(paste0("correct_edger", "_", "true_sen", "_fit.RDS"))
df$symbol <- rownames(df)
df$sig <- df$FDR <= 0.05

true_sen <- df
var <- read.delim("var_de_aucell_mg_exc.txt", header = FALSE)
var <- var[, 1:3]
var <- unique(var)

i <- 5
file_path <- paste0("correct_edger", "_", var[i, 2], "_fit.RDS")
df <- readRDS(file_path)
df$test <- var[i, 2]
df$symbol <- rownames(df)

all2 <- df
all2$sig <- all2$FDR <= 0.05
predicted_sen <- all2

all <- merge(true_sen, predicted_sen, by = c("symbol"), all.y = TRUE)

correlations <- all %>%
  summarise(
    test_result = list(cor.test(logFC.x, logFC.y, method = "spearman", use = "complete.obs"))
  ) %>%
  mutate(
    rho    = sapply(test_result, function(x) x$estimate),
    p_value = sapply(test_result, function(x) x$p.value),
    label  = sprintf("rho = %.2f, p = %.3f", rho, p_value)
  )

d <- ggplot(all, aes(x = logFC.x, y = logFC.y)) +
  geom_point(color = "grey40", alpha = 0.7) +
  theme_bw(base_size = 14) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "True Senescence LogFC",
    y = "Predicted Senescence\nLogFC"
  ) +
  theme(
    text       = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20)
  ) +
  geom_text(
    data = correlations,
    aes(label = label),
    x = 0.95 * max(all$logFC.x, na.rm = TRUE),
    y = 0.95 * max(all$logFC.y, na.rm = TRUE),
    hjust = 1, vjust = 1,
    size = 10
  )

all3 <- ((b + h + d) +
  plot_layout(axes = "collect"))

################################################################################################################
# FINAL PLOT
################################################################################################################

final <- ((n / predicted_row1 / true_row1 / all3) +
  plot_annotation(tag_levels = "A")) &
  theme(plot.tag = element_text(size = 30))

ggsave("aucell_validation_plot.pdf",
       final, width = 18, height = 15)

