library(data.table)
library(ggplot2)
library(dplyr)


senCID = readRDS("all_sc_senCID_output.RDS") #senCID scores

mg = #aucell output at chosen threshold
exc = #aucell output at chosen threshold
int = #aucell output at chosen threshold
opc = #aucell output at chosen threshold
ast = #aucell output at chosen threshold
oli = #aucell output at chosen threshold

aucell_results = rbind(mg, exc, int, nonneu, opc, ast, oli)
aucell_results$V1 = rownames(aucell_results)

final_all = merge(aucell_results, senCID, by = "V1")


##check what reSID - which chosen and then run fishers test for each cell type
library(data.table)
library(ggpubr)

dirs <- list.dirs(path = "/SenCID/", full.names = TRUE, recursive = FALSE)

dirs = dirs[-1]

df_final=c()
for (x in dirs){
	print(x)
	 # files <- list.files(dirs)
	 resSID = fread(paste0(x, "/recSID.csv"))
	 df_final = rbind(resSID,df_final)
}

aucell = aucell_results[,c("celltype", "percent")]

aucell$V1 = rownames(aucell)

all = merge(df_final, aucell, by = "V1")
all = as.data.frame(all)
all2 = merge(all[,c(1,2)], final_all, by.x = c("V1","RecSID"), by.y = c("V1", "SID"))

all3 = all2[all2$celltype != "NonNeu", ]

correlations <- all3 %>%
  summarise(
    cor_test = list(cor.test(Decision, sen_list, method = "spearman")),
    rho = cor_test[[1]]$estimate[[1]],  # Extract rho properly
    p_value = cor_test[[1]]$p.value     # Extract p-value properly
  ) %>%
  ungroup() %>%
  mutate(
    fdr = p.adjust(p_value, method = "fdr"),
    label = sprintf("rho = %.2f, p < 2.23e-308", rho)
  )


# Create the plot
g = ggplot(all3, aes(x = Decision, y = sen_list)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_label(data = correlations, 
             aes(x = -Inf, y = Inf, label = label),
             hjust = -0.1, vjust = 1.5,
             size = 4, fill = "white", color = "black") +
  theme_bw() +
  theme(
    strip.text.x = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ) +
  labs(x = "SenCID Decision", y = "Senescence Score (AUCell)") +
  ylim(-0.01, 0.1)

ggsave("aucell_sencid_sc_rec_only.pdf", g, width = 11, height = 7)
