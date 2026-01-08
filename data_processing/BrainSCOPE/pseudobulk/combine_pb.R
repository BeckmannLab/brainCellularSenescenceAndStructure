##combine pseudobulk


######################
###combine pseudobulk
#######################
library(qs)

intersectSCE <- function(...) {
  # Get all SCE objects passed to the function
  sce_list <- list(...)
  
  # Find common features (genes) across all datasets
  common_features <- Reduce(intersect, lapply(sce_list, rownames))
  
  # Subset all objects to common features
  sce_list <- lapply(sce_list, function(sce) sce[common_features,])
  
  # Combine all objects
  combined <- do.call(cbind, sce_list)
  
  # Add batch information using file names
  combined$batch <- rep(names(sce_list), sapply(sce_list, ncol))
  
  return(combined)
}

resfinal = foreach(i = percentages) %dopar% {
    print(i)
    files <- list.files(path = "../all_pseudo_bulk/", pattern = i, full.names = TRUE)   
    
    # List to store individual SCE objects
    sce_list <- list()
    
    for (x in files) {
        # Load the SCE object
        sce <- qread(x)
        
        # Check if the SCE object has exactly 12 assays
        if (length(assays(sce)) == 12) {
            mcols(sce) <- NULL
            
            # Keep only specified columns in colData
            keep_cols <- c("Cohort", "Individual_ID")
            colData(sce) <- colData(sce)[, keep_cols, drop=FALSE]
            
            # Add to the list with file name as the list name
            sce_list[[basename(x)]] <- sce
        }
    }

    # Use the intersectSCE function to combine all datasets
    if (length(sce_list) > 0) {
        final_pb <- do.call(intersectSCE, sce_list)
    } else {
        final_pb <- NULL
    }
    
    qsave(final_pb, paste0(i,".qs"))
}


#add aggr
for (i in percentages){
  print(i)
  pb = readRDS(paste0(i,".RDS"))
  obj=pb
  tmp=do.call(rbind,metadata(obj)[grepl("aggr_means$", names(metadata(obj)))])
  tmp2 = metadata(obj)$agg_pars
  metadata(obj) = list()
  metadata(obj)$aggr_means = NA
  metadata(obj)$aggr_means =tmp
  metadata(obj)$agg_pars <- tmp2
  qsave(obj, paste0("pseudobulk_with_agg",i,".qs"))
}

