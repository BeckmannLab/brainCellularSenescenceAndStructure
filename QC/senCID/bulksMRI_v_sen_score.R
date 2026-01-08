library(lme4)
library(lmerTest)
library(data.table)
str_data_summarized=readRDS("str_data_summarized_JAN172024.RDS")
str_data_summarized_weights=str_data_summarized$weights
str_data_summarized=str_data_summarized$E
cov_str_summarized=readRDS("cov_str_summarized_JAN172024.RDS")
cov_str_summarized=cov_str_summarized[match(colnames(str_data_summarized),cov_str_summarized$subject_id),]
cov_str_summarized$SurfaceHoles=str_data_summarized["SurfaceHoles",]
cov_str_summarized$rhSurfaceHoles=str_data_summarized["rhSurfaceHoles",]
cov_str_summarized$lhSurfaceHoles=str_data_summarized["lhSurfaceHoles",]


SID1 = fread("SID1.csv")
SID1$SID = "SID1"

SID2 = fread("SID2.csv")
SID2$SID = "SID2"

SID3 = fread("SID3.csv")
SID3$SID = "SID3"

SID4 = fread("SID4.csv")
SID4$SID = "SID4"

SID5 = fread("SID5.csv")
SID5$SID = "SID5"

SID6 = fread("SID6.csv")
SID6$SID = "SID6"

all = rbind(SID1, SID2, SID3, SID4, SID5, SID6)

cov=readRDS("bulk_meta_info.RDS")

df = merge(cov, all, by.x = "sample_id", by.y = "V1")


SIDS = c("SID1", "SID2", "SID3", "SID4", "SID5","SID6")
SIDS = c("SID4")

all_results = c()
for (i in unique(rownames(str_data_summarized))) {
  sub_df = str_data_summarized[rownames(str_data_summarized) == i,]
  sub_df = as.data.frame(sub_df)
  colnames(sub_df) = "image_feature"
  sub_df$subject_id = rownames(sub_df)

  for (x in SIDS) {
    sub_df2 = df[df$SID == x,]
    
    cmc_df = merge(sub_df, sub_df2, by.y = "subject_id")
    
    # Perform correlation test
   model <- lmer(SID_Score ~ image_feature + (1|subject_id) + (1|sex) + (1| disease_state) , data = cmc_df)
   model_summary = summary(model)
   p_values <- model_summary$coefficients[, 5] 
   p_values_fdr <- p.adjust(p_values, method = "fdr")

    results <- data.frame(
      feature = i,
      SIDS = x,
      p_value = p_values["image_feature"],
      fdr = p_values_fdr["image_feature"] 
      # p_value = cor_df$p.value
    )

    all_results <- rbind(all_results, results)
  }
}
sig_all_results = all_results[which(all_results$fdr <= 0.05),]
