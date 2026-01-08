library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=30)
library(R.matlab)
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)
library(rrvgo)
library(Matrix.utils)

annotate_GOterms_DE <- function(mel,sample_name,ensGene=TRUE,rrvgo=TRUE,gene.map=gene.map, useFDR=FALSE){
    # Performs GO annotation of the genes of interest
    # inputs are mel, a 2 column data.frame with number of rows corresponding
    # to your sample size, where the first column is gene name and the second
    # column is a binary vector of 0s and 1s, with 0s for genes in background
    # (not significantly differentially expressed genes for example) and 1s
    # for genes of interest (significantly differentially expressed genes for
    # example); sample_name is a string containing the name of the file that 
    # will be saved with the output of the analysis; and date, a string of the
    # date of the experiment that will also be used in the name of the output
    # file. The function has no output, output is directly saved as a table.
    #gofigure open source while revigo is not 
    options(stringsAsFactors = FALSE)
    library("R.matlab")
    library(goseq)
    library(topGO)
    library(org.Hs.eg.db)
    library(Rgraphviz)
    library(rrvgo)
    library(Matrix.utils)

    net = mel
    output.annotation.file = paste(sample_name,"_module_enrichments",".txt",sep="")
    
    #make directories 
    if(rrvgo==TRUE){
        output.annotation.rrvgo = paste(sample_name,"_","rrvgo",sep="")
        if(dir.exists(output.annotation.rrvgo)){
          system(paste("rm -rf ",output.annotation.rrvgo,sep=""))
        }
        dir.create(output.annotation.rrvgo,showWarnings=FALSE)
    }


     #get the GO IDS - gene set. which gene set belong to 
    #commented when running in a parallel loop and do it before you run if doing in parallel 
     if (is.null(gene.map)){
     if(ensGene==TRUE){
        gene.map = getgo(net[,1],'hg19','ensGene')       # extracting go IDs from gene names
     }else{
         gene.map = getgo(net[,1],'hg19','geneSymbol')       # extracting go IDs from gene names
     }
 }


    gene.map = gene.map[!is.na(names(gene.map))]        # removing non existing entries
    number_of_GO = length(unique(unlist(gene.map))) #21722
    
        #safty to make sure something not loaded in
    if (exists("res")){
       remove("res")
    }

    counter = 0
    counter = counter + 1
    flush.console()
 
    x = net[,1]         # extracting gene names
    a = rep(0,length(x))    # creating vector of length module names filled with 0s and name of row gene name

    names(a) = x
    a[net[,2] == 1] = 1     # changing indices in this vector at module positions to 1, vector used in definition of GO element

   ##BP = biological processes 
    ips_BP = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("BP"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)
    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_BP, test.stat)   # doing Fisher statistics on GO graph
    res.final_BP = GenTable(ips_BP, classic = res.fisher, topNodes=length(ips_BP@graph@nodes))   # selecting top 200 of res.fisher


    #mf = molecular functinos 
    ips_MF = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("MF"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)
    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_MF, test.stat)   # doing Fisher statistics on GO graph
    res.final_MF = GenTable(ips_MF, classic = res.fisher, topNodes=length(ips_MF@graph@nodes))   # selecting top 200 of res.fisher

    #cc = cellular components 
    ips_CC = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("CC"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)
    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_CC, test.stat)   # doing Fisher statistics on GO graph
    res.final_CC = GenTable(ips_CC, classic = res.fisher, topNodes=length(ips_CC@graph@nodes))   # selecting top 200 of res.fisher
    

    res=cbind(rbind(res.final_BP,res.final_MF,res.final_CC))

    res.mod = cbind(res, res[,"Significant"]/res[,"Expected"],p.adjust(res[,"classic"],method="BH",n=number_of_GO))  # adding a column of ratio of significant over expected (fold enrichment)
    names(res.mod)[c(7,8)] = c("fold_enrichment","BH")
    res.mod = res.mod[,c( "GO.ID", "Term", "Annotated", "Significant", "Expected", "fold_enrichment", "classic","BH")]

    write.table(res.mod, output.annotation.file, sep="\t", quote=FALSE, row.names=FALSE)

    if (useFDR==TRUE){var="BH"}else{var="classic"}
    
    if(rrvgo==TRUE){
    res.mod[is.na(res.mod[,var]),var] = 10^-300 
    tmp = res.mod[res.mod[,var] <=0.05,c('GO.ID',var)]
    go_analysis <- tmp
    simMatrix_BP <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
    simMatrix_MF <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")
    simMatrix_CC <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")

    combo_simMat <- rBind.fill(simMatrix_BP, simMatrix_CC)
    combo_all <- rBind.fill(combo_simMat, simMatrix_MF)
    combo_all <- na.omit(combo_all)
    go_analysis[,var] <- as.numeric(go_analysis[,var])
    scores <- setNames(-log10(go_analysis[,var]), go_analysis$GO.ID)
    reducedTerms <- reduceSimMatrix(combo_all,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

	# plot_place = "plots/"
    pdf(paste(plot_directory,output.annotation.rrvgo, "treemap_output.annotation.pdf", sep=""))
    scatterPlot(combo_all,reducedTerms)
    wordcloudPlot(reducedTerms,min.freq=1, colors="black" )
    heatmapPlot(combo_all,reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=6 )
    treemapPlot(reducedTerms)
    dev.off()
	}
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


setwd("path/")
var = structure(list(V1 = c("lbp_06.21.24_pseudobulk_1percent.RDS", 
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS", 
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS", 
"lbp_06.21.24_pseudobulk_1percent.RDS", "lbp_06.21.24_pseudobulk_1percent.RDS", 
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS", 
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS", 
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_5percent.RDS", 
"lbp_06.21.24_pseudobulk_5percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS", 
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS", 
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS", 
"lbp_06.21.24_pseudobulk_20percent.RDS", "lbp_06.21.24_pseudobulk_20percent.RDS", 
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS", 
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS", 
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_30percent.RDS", 
"lbp_06.21.24_pseudobulk_30percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS", 
"lbp_06.21.24_pseudobulk_10percent.RDS", "lbp_06.21.24_pseudobulk_10percent.RDS"
), V2 = c("1percent", "1percent", "1percent", "1percent", "1percent", 
"1percent", "1percent", "5percent", "5percent", "5percent", "5percent", 
"5percent", "5percent", "5percent", "20percent", "20percent", 
"20percent", "20percent", "20percent", "20percent", "20percent", 
"30percent", "30percent", "30percent", "30percent", "30percent", 
"30percent", "30percent", "10percent", "10percent", "10percent", 
"10percent", "10percent", "10percent", "10percent"), V3 = c("pb_1percent", 
"pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", "pb_1percent", 
"pb_1percent", "pb_5percent", "pb_5percent", "pb_5percent", "pb_5percent", 
"pb_5percent", "pb_5percent", "pb_5percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_20percent", "pb_20percent", 
"pb_20percent", "pb_20percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_30percent", "pb_30percent", "pb_30percent", 
"pb_30percent", "pb_10percent", "pb_10percent", "pb_10percent", 
"pb_10percent", "pb_10percent", "pb_10percent", "pb_10percent"
), V4 = c("Ast_TRUE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", 
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE", "Ast_FALSE", 
"Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", 
"OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", "MG_TRUE", "NonNeu_TRUE", 
"Oli_TRUE", "OPC_TRUE", "Ast_FALSE", "Exc_TRUE", "Int_TRUE", 
"MG_TRUE", "NonNeu_TRUE", "Oli_TRUE", "OPC_TRUE"), V5 = c("Ast_FALSE", 
"Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", 
"OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", 
"NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", 
"Int_FALSE", "MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE", 
"Ast_TRUE", "Exc_FALSE", "Int_FALSE", "MG_FALSE", "NonNeu_FALSE", 
"Oli_FALSE", "OPC_FALSE", "Ast_TRUE", "Exc_FALSE", "Int_FALSE", 
"MG_FALSE", "NonNeu_FALSE", "Oli_FALSE", "OPC_FALSE"), V6 = c("ast", 
"exc", "int", "mg", "noneu", "oli", "opc", "ast", "exc", "int", 
"mg", "noneu", "oli", "opc", "ast", "exc", "int", "mg", "noneu", 
"oli", "opc", "ast", "exc", "int", "mg", "noneu", "oli", "opc", 
"ast", "exc", "int", "mg", "noneu", "oli", "opc"), V7 = c(0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)), row.names = c(NA, -35L
), class = "data.frame")


all2 = c()
for (i in 1:nrow(var)){
    df = readRDS(paste0("path/",var[i,6],"_",var[i,2], "_fit.RDS"))
    df$celltype = var[i,6]
    df$test = var[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}

resFinal = all2
colnames(resFinal)[5] = "adj.P.Val"


sample_name = "GO_sen"
plot_directory = "path"
resFinal$DE = "notSignif"
resFinal$DE = ifelse(resFinal$adj.P.Val <= 0.05, "Signif", "notSignif")
cell_type = unique(resFinal$celltype)

percent = c("30percent", "10percent", "1percent", "30percent", "5percent", "5percent", "1percent")
var = cbind(cell_type, percent)
var = as.data.frame(var)


foreach(i = 1:nrow(var)) %do% {
    df = resFinal[which(resFinal$celltype == var[i,1] & resFinal$test == var[i,2]),]
    df = df[df$logFC >= 0,]
    for_gsea_df = format_for_gsea(data.frame(gene = df$symbol, sig = df$DE))
    for_gsea_df$DE_new = ifelse(for_gsea_df$DE == "Signif", 1, 0)
    for_gsea_df = for_gsea_df[,c("gene_ID", "DE_new")]
    sample_name = paste0("GO_sen_results","_",var[i,1],"_", var[i,2],"_", "positive")
    gene.map = getgo(for_gsea_df$gene_ID,'hg19','geneSymbol')
    print(paste(i))
    if (any(for_gsea_df$DE_new == 1)) {
    try(annotate_GOterms_DE(mel=for_gsea_df,sample_name=sample_name,gene.map=gene.map,rrvgo=TRUE,useFDR=TRUE,ensGene==FALSE))
}
else {
	print("no sig")
}	
}


foreach(i = 1:nrow(var)) %do% {
    df = resFinal[which(resFinal$celltype == var[i,1] & resFinal$test == var[i,2]),]
    df = df[df$logFC < 0,]
    for_gsea_df = format_for_gsea(data.frame(gene = df$symbol, sig = df$DE))
    for_gsea_df$DE_new = ifelse(for_gsea_df$DE == "Signif", 1, 0)
    for_gsea_df = for_gsea_df[,c("gene_ID", "DE_new")]
    sample_name = paste0("GO_sen_results","_",var[i,1],"_", var[i,2],"_", "negative")
    gene.map = getgo(for_gsea_df$gene_ID,'hg19','geneSymbol')
    print(paste(i))
    if (any(for_gsea_df$DE_new == 1)) {
    try(annotate_GOterms_DE(mel=for_gsea_df,sample_name=sample_name,gene.map=gene.map,rrvgo=TRUE,useFDR=TRUE,ensGene==FALSE))
}
else {
    print("no sig")
}   
}


#####
#pull results

library(data.table)
library(stringr)
flist <- Sys.glob("path/*txt")
myres <- c()
failures <- c()
for (i in flist){
  print(i)
  data <- fread(i)
  celltype = gsub("path/", "", i)
  split_string <- strsplit(celltype, "_")
  celltype2 <- split_string[[1]][1]
  direction <- split_string[[1]][3]

  data$direction = direction
  data$celltype = celltype2
  myres = rbind(data,myres )
}


saveRDS(myres,"sen_sig_GO_by_features_7.8.24.RDS" )

