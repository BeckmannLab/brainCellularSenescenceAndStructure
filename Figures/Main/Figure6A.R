########
#figure6A
##########
ml R/4.4.1
R

library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(scales)

exc_nd = #aucell results
int_nd = #aucell results
opc_nd = #aucell results
vas_nd = #aucell results
ast_nd = #aucell results
oli_nd = #aucell results
micro_nd = #aucell results

meta= #meta data

cov2 = meta[,c("batch","RL#", "age", "numerical_age")]
cov2_nd = unique(cov2)

V1 = c("exc_nd", "int_nd", "opc_nd", "oli_nd", "ast_nd", "vas_nd", "micro_nd")
V3 = c("exc", "int", "opc", "oli", "ast", "nonneu", "mg")
var = cbind(V1, V3)

var = as.data.frame(var)

exc = list()
int = list()
opc = list()
oli = list()
ast = list()
nonneu = list()
mg = list()
all = c()

for (i in 1:nrow(var)) {
        #get ids
        nd = get(var[i,1])

        ##format subject IDs
        nd <- nd %>%
        mutate(Sample_ID = str_extract(cell, "(?<=-).*"))

        #percentage
        df_nd <- nd %>%
        group_by(Sample_ID) %>%
        summarize(
        percentage_sen = mean(sen == TRUE) * 100
        )

        ##add covariates
        all_nd = merge(df_nd, cov2_nd, by.x = "Sample_ID", by.y = "batch")
  
        #format
        colnames(all_nd)[5] = c("Age")

        #subset
        all_nd_subset = all_nd[,c("Sample_ID", "percentage_sen", "Age")]
        all_nd_subset$ds = "nd"

        #rbind it
        print(i)
        myres= all_nd_subset
        myres$celltype = var[i, 2]

        all = rbind(myres, all)
        assign(var[i, 2],myres)
}

all_mg = all[which(all$celltype == "mg"),]
all_exc = all[which(all$celltype == "exc"),]
all_nonneu = all[which(all$celltype == "nonneu"),]
all_ast = all[which(all$celltype == "ast"),]
all_oli = all[which(all$celltype == "oli"),]
all_opc = all[which(all$celltype == "opc"),]
all_int = all[which(all$celltype == "int"),]

all$celltype = gsub("mg", "Microglia", all$celltype)
all$celltype = gsub("exc", "Excitatory Neurons", all$celltype)
all$celltype = gsub("ast", "Astrocytes", all$celltype)
all$celltype = gsub("oli", "Oligodendrocytes", all$celltype)
all$celltype = gsub("opc", "Oligodendrocyte progenitor cells", all$celltype)
all$celltype = gsub("int", "Inhibitory Neurons", all$celltype)

all <- all[all$celltype != "nonneu", ]

#plots
b = ggplot(all, aes(x = Age, y = percentage_sen)) +
  geom_point(size = 2.5, alpha = 0.5) +  # Add points with size 2.5
  scale_x_continuous(name = "Age", 
                     limits = c(-5, max(all$Age, na.rm = TRUE))) +  # Remove padding
  scale_y_continuous(name = "Percent Senescent", labels = percent_format(scale = 1)) +
  geom_smooth(method = loess) + 
  theme_bw() +  # Use a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 22),  # Increase x-axis title size
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 18),
    legend.position = "none"
  ) + 
  facet_wrap(~celltype, scales = "free")

ggsave("nd_celltype_cor_age_sens.pdf", b, height = 10, width = 13)
