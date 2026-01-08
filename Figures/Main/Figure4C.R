### og
library(qvalue)
library(data.table)
library(seriation)
library(corrplot)
library(ggplot2)

## original dataset: liv aucell
var  = read.delim("var_de_aucell.txt", header = FALSE)
var2 = read.delim("var_de_aucell_10.txt", header = FALSE)
var  = rbind(var, var2)
var_for_og = var[c(1, 9, 10, 25, 34, 28), ]

all2 = c()
for (i in 1:nrow(var_for_og)) {
    df = readRDS(paste0(
        "correct_edger_",
        var_for_og[i, 6], "_", var_for_og[i, 2], "_fit.RDS"
    ))
    df$celltype = var_for_og[i, 6]
    df$test     = var_for_og[i, 2]
    df$symbol   = rownames(df)
    all2        = rbind(df, all2)
}

liv_aucell = all2

## psychencode (brainSCOPE)
var = readRDS("var_file.RDS")

var_psychencode = var[c(30, 25, 27, 14, 16, 29), ]

all2 = c()
for (i in 1:nrow(var_psychencode)) {
    df = readRDS(paste0(
        "correct_edger_",
        var_psychencode[i, 6], "_", var_psychencode[i, 2], "_fit.RDS"
    ))
    df$celltype = var_psychencode[i, 6]
    df$test     = var_psychencode[i, 2]
    df$symbol   = rownames(df)
    all2        = rbind(df, all2)
}

psychencode_aucell = all2

## pm dataset
var  = read.delim("var_de_aucell.txt", header = FALSE)
var2 = read.delim("var_de_aucell_10.txt", header = FALSE)
var  = rbind(var, var2)
var_pm = var[c(13, 2, 10, 4, 29, 21), ]

all2 = c()
for (i in 1:nrow(var_pm)) {
    df = readRDS(paste0(
        "correct_edger_",
        var_pm[i, 6], "_", var_pm[i, 2], "_fit.RDS"
    ))
    df$celltype = var_pm[i, 6]
    df$test     = var_pm[i, 2]
    df$symbol   = rownames(df)
    all2        = rbind(df, all2)
}
pm_aucell = all2

## correlations

# og to psychencode
var_for_og2      = var_for_og[, c(2, 6)]
var_psychencode2 = var_psychencode[, c(2, 6)]

var_for_og2      = var_for_og2[order(var_for_og2$V6), ]
var_psychencode2 = var_psychencode2[order(var_psychencode2$v6), ]

all_var = cbind(var_psychencode2, var_for_og2)

myres = c()

for (x in 1:nrow(all_var)) {
    liv_aucell_subset = liv_aucell[
        which(liv_aucell$celltype == all_var[x, 4] & liv_aucell$test == all_var[x, 3]),
    ]
    psychencode_aucell_subset = psychencode_aucell[
        which(psychencode_aucell$celltype == all_var[x, 2] & psychencode_aucell$test == all_var[x, 1]),
    ]
    mer  = merge(liv_aucell_subset, psychencode_aucell_subset, by = "symbol")
    rho  = cor.test(mer$logFC.x, mer$logFC.y, method = "spearman")
    p_val = rho$p.value
    fdr   = p.adjust(p_val, method = "fdr")
    add   = data.table(
        celltype  = all_var[x, 2],
        spearmans = rho$estimate,
        p         = p_val,
        fdr       = fdr
    )
    myres = rbind(myres, add)
    print(add, row.names = FALSE, col.names = "none")
}

og_psychencode       = myres
og_psychencode$test  = "og_psychencode"

# og to pm 
var_for_og2 = var_for_og[, c(2, 6)]
var_pm2     = var_pm[, c(2, 6)]
colnames(var_pm2) = c("v2", "v6")

var_for_og2 = var_for_og2[order(var_for_og2$V6), ]
var_pm2     = var_pm2[order(var_pm2$v6), ]

all_var = cbind(var_pm2, var_for_og2)

myres = c()

for (x in 1:nrow(all_var)) {
    liv_aucell_subset = liv_aucell[
        which(liv_aucell$celltype == all_var[x, 4] & liv_aucell$test == all_var[x, 3]),
    ]
    pm_aucell_subset = pm_aucell[
        which(pm_aucell$celltype == all_var[x, 2] & pm_aucell$test == all_var[x, 1]),
    ]
    mer  = merge(liv_aucell_subset, pm_aucell_subset, by = "symbol")
    rho  = cor.test(mer$logFC.x, mer$logFC.y, method = "spearman")
    p_val = rho$p.value
    fdr   = p.adjust(p_val, method = "fdr")
    add   = data.table(
        celltype  = all_var[x, 2],
        spearmans = rho$estimate,
        p         = p_val,
        fdr       = fdr
    )
    myres = rbind(myres, add)
    print(add, row.names = FALSE, col.names = "none")
}

og_pm      = myres
og_pm$test = "og_pm"

### plot

all_cors = rbind(og_pm, og_psychencode)

all_cors$new_celltype = all_cors$celltype
all_cors$new_celltype = gsub("Astro", "Astrocytes",                    all_cors$new_celltype)
all_cors$new_celltype = gsub("ast",   "Astrocytes",                    all_cors$new_celltype)
all_cors$new_celltype = gsub("exc",   "Excitatory Neurons",            all_cors$new_celltype)
all_cors$new_celltype = gsub("Exc",   "Excitatory Neurons",            all_cors$new_celltype)
all_cors$new_celltype = gsub("int",   "Inhibitory Neurons",            all_cors$new_celltype)
all_cors$new_celltype = gsub("Int",   "Inhibitory Neurons",            all_cors$new_celltype)
all_cors$new_celltype = gsub("micro", "Microglia",                     all_cors$new_celltype)
all_cors$new_celltype = gsub("Micro", "Microglia",                     all_cors$new_celltype)
all_cors$new_celltype = gsub("Oligo", "Oligodendrocytes",              all_cors$new_celltype)
all_cors$new_celltype = gsub("oli",   "Oligodendrocytes",              all_cors$new_celltype)
all_cors$new_celltype = gsub("OPC",   "Oligodendrocyte\nProgenitor Cells", all_cors$new_celltype)
all_cors$new_celltype = gsub("opc",   "Oligodendrocyte\nProgenitor Cells", all_cors$new_celltype)

all_cors$new_celltype = gsub("Excitatory Neuronsitatory Neurons", "Excitatory Neurons", all_cors$new_celltype)
all_cors$new_celltype = gsub("Microgliaglia",                     "Microglia",          all_cors$new_celltype)
all_cors$new_celltype = gsub("mg",                                "Microglia",          all_cors$new_celltype)

all_cors$test = gsub("og_pm",          "Brain Bank matched samples", all_cors$test)
all_cors$test = gsub("og_psychencode", "brainSCOPE",                all_cors$test)

g = ggplot(all_cors,
           aes(x = reorder(new_celltype, -spearmans),
               y = spearmans,
               color = test)) +
  geom_point(size = 6) +
  theme_bw() +
  labs(x     = "Cell Type",
       y     = "Correlation (Spearman's \u03C1)",
       color = "Dataset") +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  theme(
    legend.position       = c(0.95, 0.95),
    legend.justification  = c("right", "top"),
    legend.box.background = element_rect(colour = "black"),
    axis.title.x          = element_text(size = 19),
    axis.title.y          = element_text(size = 19),
    axis.text.x           = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y           = element_text(size = 18),
    legend.text           = element_text(size = 14),
    legend.title          = element_text(size = 17)
  )

ggsave("all_liv_sen_correlation.pdf",
       g, width = 7, height = 7)
