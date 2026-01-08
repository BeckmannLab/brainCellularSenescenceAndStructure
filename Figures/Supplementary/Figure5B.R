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

###plot

df$celltype = factor(df$celltype)

df2 <- df %>%
  group_by(celltype) %>%
  wilcox_test(percentage_sen ~ mymet_PD, alternative = "two.sided", detailed = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


df3 = df[which(df$celltype!= "nonneu"),]

df3$celltype = gsub("ast", "Ast", df3$celltype)
df3$celltype = gsub("ast", "Ast", df3$celltype)
df3$celltype = gsub("exc", "Exc", df3$celltype)
df3$celltype = gsub("int", "Inh", df3$celltype)
df3$celltype = gsub("mg", "MG", df3$celltype)
df3$celltype = gsub("oli", "Oli", df3$celltype)
df3$celltype = gsub("opc", "OPC", df3$celltype)


df3$celltype = factor(df3$celltype, levels = c("Ast", "Exc", "Inh", "MG", "Oli", "OPC"))

d= ggboxplot(df3, x = "celltype", y = "percentage_sen",
          color = "mymet_PD", palette = "jco") +
  scale_x_discrete(name = "Cell Type") +
  scale_y_continuous(name = "Percent Senescence", 
                     labels = scales::percent_format(scale = 1),
                     limits = c(0, 100)) +
  labs(color = "PD") +
  stat_summary(fun = max, geom = "text", aes(label = "ns"), 
               vjust = -0.5, color = "black") +
  theme_pubr() +
  theme(axis.text.x = element_text(hjust = 1))

ggsave("percentage_sen_individual_boxplot_PD_with_fdr.pdf",d)
