ml R/4.4.1
R

library(data.table)
library(SingleCellExperiment)
library(zellkonverter)
library(spam)
library(spam64)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(Seurat)
library(AUCell)
library(dplyr)
library(Matrix)
library(harmony)
library(BiocParallel)
Sys.setenv(OMP_NUM_THREADS = 48)
library(readxl)
library(GSEABase)
library(HDF5Array)
library(dreamlet)
library(foreach)
library(stringr)
library(tidyr)
library(qs)

#psychgen
PEC2_sample=fread("PEC2_sample_metadata.txt",data.table=FALSE)
PEC2_mapping=as.data.frame(readxl::read_excel("PEC2_sample_mapping.xlsx"))

#fix for format
PEC2_sample[which(PEC2_sample$Cohort == "Girgenti-snMultiome"),"Cohort"] = "Girgenti-multiome"
PEC2_mapping[which(PEC2_mapping$Cohort == "Girgenti-snMultiome"),"Cohort"] = "Girgenti-multiome"
PEC2_sample[which(PEC2_sample$Individual_ID == "HSB340"),"Individual_ID"] = "HSB340_DFC"
duplicates = c("HSB189", "HSB340", "UMB5297", "RT00389N")

##all meta
snRNAseqPFC_BA10_biospecimen_metadata <- fread("snRNAseqPFC_BA10_biospecimen_metadata.csv")
clinical = fread("ROSMAP_clinical.csv")
ROSMAP_biospecimen_metadata = fread("ROSMAP_biospecimen_metadata.csv")
ROSMAP_assay_scrnaSeq_metadata = fread("ROSMAP_assay_scrnaSeq_metadata.csv")
filtered.colMetadata <- read.delim("filtered_column_metadata.txt")

##format
meta = merge(snRNAseqPFC_BA10_biospecimen_metadata, clinical, by = "projid", all.x = TRUE)
folder = "ROSMAP"
meta$PMI <- round(meta$pmi, digits = 1)
meta$Biological_Sex = ifelse(meta$msex == 0, "female", "male")
meta$Age_death = gsub("\\+", "", meta$age_death)
meta$Age = gsub("\\+", "", meta$age_death)
meta$Age = round(as.numeric(meta$Age), digits = 0)
meta$Disorder = case_when(
  meta$dcfdx_lv == 1 ~ "control",
  meta$dcfdx_lv %in% 2:3 ~ "cognitive impairment",
  meta$dcfdx_lv %in% 4:6 ~ "Alzheimers/dementia",
  TRUE ~ NA_character_
)
meta = as.data.frame(meta)
meta$Cohort = "ROSMAP"
meta$Individual_ID = meta$projid
meta$Ancestry = meta$race
meta$"1000G_ancestry" = NA
meta$Genotype_data = "RNA" 
meta$"snATAC-Seq" = NA
meta$pH = NA
meta$RIN =  NA 
meta$Notes = NA

final_meta = meta[,c("Cohort", "Individual_ID", "Biological_Sex", "Age_death", "Disorder", "1000G_ancestry", "Genotype_data", "snATAC-Seq", "PMI", "pH", "RIN", "Notes", "Age")]
meta = merge(final_meta, filtered.colMetadata, by.x = "Individual_ID", by.y = "projid")
meta$cell_id <- paste0(gsub("\\.\\d+$", "", meta$TAG), "_", meta$Individual_ID)
rownames(meta) = meta$cell_id

#expression
filtered.counts <- readMM("filtered_count_matrix.mtx")
rownames(filtered.counts) <- readLines("filtered_gene_row_names.txt")
filtered.colMetadata$cell_id <- paste0(gsub("\\.\\d+$", "", filtered.colMetadata$TAG), "_", filtered.colMetadata$projid)
colnames(filtered.counts) = filtered.colMetadata$cell_id

##format
countMat = filtered.counts
countMat <- as(countMat, "dgCMatrix")

order <- match(rownames(meta), colnames(countMat))
countMat<- countMat[, order]

#convert to sce
sce <- SingleCellExperiment(assays = list(counts=countMat), colData=meta)

colData(sce)$new_celltype <- ifelse(colData(sce)$broad.cell.type %in% 
                                      c("Ex"), "Exc", 
                                      ifelse(colData(sce)$broad.cell.type %in% 
                                             c("In"), "Int", 
                                             ifelse(colData(sce)$broad.cell.type %in% 
                                             c("Opc"), "OPC",
                                             ifelse(colData(sce)$broad.cell.type %in% 
                                             c("Oli"), "Oligo",
                                             ifelse(colData(sce)$broad.cell.type %in% 
                                             c("Ast"), "Astro",
                                             ifelse(colData(sce)$broad.cell.type %in% 
                                             c("Mic"), "Micro",
                                             colData(sce)$broad.cell.type))))))

saveRDS(sce,"sce_all_rosmap.RDS")
