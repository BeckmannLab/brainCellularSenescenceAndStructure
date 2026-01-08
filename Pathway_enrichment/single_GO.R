library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
library(foreach)
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(rrvgo)

##functions
annotate_GOterms_DE <- function(
	  mel,
	  sample_name,
	  ensGene = TRUE,
	  rrvgo = TRUE,
	  gene.map = NULL,          # <- default should be NULL
	  useFDR = FALSE,
	  plot_directory = ".",
	  out_dir = "."
	){
	  options(stringsAsFactors = FALSE)
	  library(goseq)
	  library(topGO)
	  library(org.Hs.eg.db)
	  library(Rgraphviz)
	  library(rrvgo)
	  library(Matrix.utils)

	  net <- mel

	  # --- ensure dirs exist ---
	  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	  dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)

	  # table path goes to out_dir
	  output.annotation.file <- file.path(out_dir, paste0(sample_name, "_module_enrichments.txt"))

	  # figures go under plot_directory/<sample_name>_rrvgo
	  if (rrvgo) {
	    plots_dir <- file.path(plot_directory, paste0(sample_name, "_rrvgo"))
	    if (dir.exists(plots_dir)) unlink(plots_dir, recursive = TRUE, force = TRUE)
	    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
	  }

	  # --- GO mapping if not provided ---
	  if (is.null(gene.map)) {
	    if (ensGene) {
	      gene.map <- getgo(net[,1], 'hg19', 'ensGene')
	    } else {
	      gene.map <- getgo(net[,1], 'hg19', 'geneSymbol')
	    }
	  }
	  gene.map <- gene.map[!is.na(names(gene.map))]

	  number_of_GO <- length(unique(unlist(gene.map)))

	  x <- net[,1]
	  a <- rep(0, length(x))
	  names(a) <- x
	  a[net[,2] == 1] <- 1
	  allGenes <- as.factor(a)

	  # --- topGO for BP/MF/CC ---
	  mk_tgd <- function(ont) new("topGOdata",
	                              description = "GO enrichment",
	                              ontology = ont,
	                              allGenes = allGenes,
	                              geneSel = function(g) g == 1,
	                              nodeSize = 10,
	                              annot = annFUN.gene2GO,
	                              gene2GO = gene.map)

	  ips_BP <- mk_tgd("BP")
	  ips_MF <- mk_tgd("MF")
	  ips_CC <- mk_tgd("CC")

	  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

	  res_BP <- GenTable(ips_BP, classic = getSigGroups(ips_BP, test.stat),
	                     topNodes = length(ips_BP@graph@nodes))
	  res_MF <- GenTable(ips_MF, classic = getSigGroups(ips_MF, test.stat),
	                     topNodes = length(ips_MF@graph@nodes))
	  res_CC <- GenTable(ips_CC, classic = getSigGroups(ips_CC, test.stat),
	                     topNodes = length(ips_CC@graph@nodes))

	  res <- rbind(res_BP, res_MF, res_CC)
	  res.mod <- cbind(
	    res,
	    res[,"Significant"] / res[,"Expected"],
	    p.adjust(res[,"classic"], method = "BH", n = number_of_GO)
	  )
	  names(res.mod)[c(7,8)] <- c("fold_enrichment","BH")
	  res.mod <- res.mod[, c("GO.ID","Term","Annotated","Significant","Expected","fold_enrichment","classic","BH")]

	  # --- write the table to out_dir ---
	  write.table(res.mod, output.annotation.file, sep = "\t", quote = FALSE, row.names = FALSE)

	  # --- rrvgo figures to plots_dir ---
	  if (rrvgo) {
	    var <- if (useFDR) "BH" else "classic"
	    # coerce p-values safely and filter
	    pnum <- suppressWarnings(as.numeric(res.mod[[var]]))
	    pnum[is.na(pnum)] <- 1e-300
	    keep <- which(pnum <= 0.05)
	    if (length(unique(res.mod$GO.ID[keep])) >= 2L) {
	      go_ids <- res.mod$GO.ID[keep]
	      scores <- stats::setNames(-log10(pnum[keep]), go_ids)

	      # one PDF with all plots in the sample's plots folder
	      pdf(file.path(plots_dir, paste0(sample_name, "_rrvgo_plots.pdf")))

	      # do each ontology separately (saves/plots into the same PDF)
	      for (ont in c("BP","MF","CC")) {
	        sim <- suppressMessages(
	          rrvgo::calculateSimMatrix(go_ids, orgdb = org.Hs.eg.db::org.Hs.eg.db, ont = ont, method = "Rel")
	        )
	        if (is.matrix(sim) && nrow(sim) >= 2 && ncol(sim) >= 2) {
	          red <- suppressMessages(
	            rrvgo::reduceSimMatrix(sim, scores, threshold = 0.7, orgdb = org.Hs.eg.db::org.Hs.eg.db)
	          )
	          rrvgo::scatterPlot(sim, red)
	          rrvgo::wordcloudPlot(red, min.freq = 1)
	          rrvgo::heatmapPlot(sim, red, annotateParent = TRUE, annotationLabel = "parentTerm", fontsize = 6)
	          rrvgo::treemapPlot(red)
	        }
	      }
	      invisible(dev.off())
	    }
	  }

	  invisible(list(table = output.annotation.file,
	                 plots_dir = if (rrvgo) plots_dir else NULL))
	}


resFinal <- #supplementary table 1 sn sMRI DE

features  <- unique(resFinal$feature)
cell_type <- unique(resFinal$assay)

new_df <- data.table()
foreach (x = features) %do% {
	print(x)
  foreach (i = cell_type) %do% {
  	print(i)
    df  <- resFinal[which(resFinal$feature == x & resFinal$assay == i), ]
    pi1 <- 1 - qvalue(df$P.Value)$pi0
    add <- data.table(feature = x, assay = i, pi1 = pi1)
    new_df <- rbind(add, new_df)
  }
}

merged_data <- merge(resFinal, new_df, by = c("feature", "assay"), all.x = TRUE)

##pi1 
length_df <- data.table()
foreach (x = features) %do% {
	print(x)
  foreach (i = cell_type) %do% {
  	print(i)
    total <- sum(merged_data$feature == x & merged_data$assay == i)
    add   <- data.table(feature = x, assay = i, total = total)
    length_df <- rbind(add, length_df)
  }
}

merged_data2 <- merge(merged_data, length_df, by = c("feature", "assay"), all.x = TRUE)
merged_data2$pi1_new <- merged_data2$total * merged_data2$pi1


##run GO 
out_dir <- "path"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

foreach (x = features) %do% {
	foreach (i = cell_type) %do% {
	df = merged_data2[which(merged_data2$feature == x & merged_data2$assay == i),]
	cutoff = unique(df$pi1_new)
	new_df = data.frame(gene_name = df$ID, DE = 0)
	df=df[order(df$P.Value, decreasing = FALSE),]
	new_df$DE[1:floor(cutoff)]=1
	sample_name = paste0("GO_bulk_results","_",x,"_", i)
	gene.map = getgo(new_df$gene_name,'hg19','geneSymbol')
	print(paste(x,"_", i))
		if (any(new_df$DE == 1)) {
		try(annotate_GOterms_DE(mel=new_df,sample_name=paste0(sample_name,"2"),gene.map=gene.map,plot_directory = "path",out_dir = out_dir, rrvgo=TRUE,useFDR=TRUE,ensGene=TRUE))
		}else {
		print("no sig")
			}	
		}
}

##pull together 
flist <- Sys.glob(file.path(out_dir, "*.txt"))

myres2 <- rbindlist(lapply(flist, function(i) {
  data <- fread(i)
  data$feature  <- gsub(paste0(out_dir, "/GO_bulk_results_scale"), "", i)
  data$feature  <- gsub("_module_enrichments.txt", "", data$feature)
  data$celltype <- str_extract(data$feature, "(?<=_)[^_]+$")
  data$feature <- sub("^scale_?", "", data$feature)
  data
}), use.names = TRUE, fill = TRUE)
