library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
library(foreach)
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(rrvgo)

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

format_for_gsea = function(df){
    #######input for this is a two column dataframe where the first column are gene symbols and second column DE status which is either Signif or notSignif. if you want to do gene IDs you need to switch out. output is two column dataframe where first column is gene symbols second column is signif or not signif
    library(assertthat)
    new_df = data.frame(gene_symbol= unique(df$gene), sig = "notSignif") ##select the columns of interest and only keep unique genes
    unique_DEGs = unique(df$gene[which(df$sig == "Signif")])
    assert_that(identical(new_df$gene[match(unique_DEGs,new_df$gene)],unique_DEGs)) #make sure they match
    
    new_df$sig[match(unique_DEGs,new_df$gene)] = "Signif" #ones that are unique become unique in the df

    assert_that(identical(sort(new_df$gene[new_df$sig=="Signif"]),sort(unique_DEGs)))
    colnames(new_df) = c("gene_ID","DE")

    ##remove empty rows 
    new_df[which(new_df$gene_ID == ""),] = NA #remove symbols that dont match
    new_df =na.omit(new_df)
    new_df
}

resFinal= #supplementary table 2 bulk sMRI DE

sample_name = "GO_bulk_results"
out_dir <- "path"


resFinal$DE = "notSignif"
resFinal$DE = ifelse(resFinal$adj.P.Val <= 0.05, "Signif", "notSignif")
features = unique(resFinal$feature)

foreach (x = features) %do% {
	df = resFinal[which(resFinal$feature == x),]
	for_gsea_df = format_for_gsea(data.frame(gene = df$gene, sig = df$DE))
	for_gsea_df$DE_new = ifelse(for_gsea_df$DE == "Signif", 1, 0)
	for_gsea_df = for_gsea_df[,c("gene_ID", "DE_new")]
	for_gsea_df$gene_id_updated=unlist(lapply(strsplit(as.character(for_gsea_df$gene_ID),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))
	for_gsea_df2 = for_gsea_df[,c(3,2)]
	sample_name = paste0("GO_bulk_results","_",x)
	gene.map = getgo(for_gsea_df2$gene_id_updated,'hg19','ensGene')
	print(x)
	if (any(for_gsea_df2$DE_new == 1)) {
	try(annotate_GOterms_DE(mel=for_gsea_df2,sample_name=sample_name,gene.map=gene.map,plot_directory = "path",out_dir = out_dir,rrvgo=TRUE,useFDR=TRUE,ensGene==TRUE))
	}
	else {
		print("no sig")
	}	
}


###
library(data.table)
flist <- Sys.glob("path/*txt")
myres <- c()
failures <- c()
for (i in flist){
  print(i)
  data <- fread(i)
  data$feature = gsub("path/", "", i)
  data$feature = gsub("_module_enrichments.txt", "", data$feature)
  myres = rbind(data,myres )
}

sig_myres = myres[which(myres$BH <= 0.05),]


