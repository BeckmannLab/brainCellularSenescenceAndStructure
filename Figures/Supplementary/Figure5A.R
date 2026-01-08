library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)
library(rstatix)

mg = #aucell results
exc = #aucell results 
int = #aucell results
opc = #aucell results
ast = #aucell results
oli = #aucell results

cells = c("mg", "exc", "int", "nonneu", "opc", "ast", "oli")

all_res= c()
for (x in cells){
  res = get(x)
  res$ID <- sub("_.*", "", rownames(res))
  df <- res %>%
  group_by(ID) %>%
  summarize(
    percentage_sen = mean(sen == TRUE) * 100
  )

  df = as.data.frame(df)
  df$celltype = x

  all_res = rbind(df, all_res)
}

cov = #metadata

cov2 = cov[,c("ids","Age", "mymet_PD","IndvID")]
colnames(cov2)[1] = "ID"
df = merge(cov2, all_res, by = "ID")

df$mymet_PD = ifelse(df$mymet_PD == "PD", TRUE,FALSE)

myres = c()
for (x in unique(df$celltype)){
subset_df = df[which(df$celltype == x),]
all_df = cor.test(subset_df$percentage_sen, subset_df$Age, alternative = "greater")
p = all_df$p.value
rho = all_df$estimate
all_f = data.table(celltype = x, p_value = p, rho = rho)
myres <- rbind(myres, all_f)
}


###plot

myres$celltype = factor(myres$celltype, levels = c("nonneu", "exc", "mg", "opc", "int", "ast", "oli"))

myres$fdr <- p.adjust(myres$p_value, method = "fdr")

myres = myres[which(myres$celltype != "nonneu"),]
myres$celltype = gsub("int", "Inh", myres$celltype)
myres$celltype = gsub("mg", "MG", myres$celltype)
myres$celltype = gsub("exc", "Exc", myres$celltype)
myres$celltype = gsub("ast", "Ast", myres$celltype)
myres$celltype = gsub("oli", "Oli", myres$celltype)
myres$celltype = gsub("opc", "OPC", myres$celltype)

myres$celltype = factor(myres$celltype, levels = c("Exc", "MG", "OPC", "Inh", "Ast", "Oli"))

d =ggplot(myres, aes(x = celltype, y = rho)) +
  geom_point(size = 3) + 
  geom_text(aes(label = round(fdr, 3)), vjust = -1, size = 4) +  # increased from 3.5 to 5
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw(base_size = 16) +  # sets a larger base text size
  labs(
    x = "Cell type",
    y = expression("Spearman's rho " * rho * " (p)")
  ) + 
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave("percentage_sen_individual_boxplot_age_with_fdr.pdf",d)
