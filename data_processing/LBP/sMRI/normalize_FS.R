
ml R/4.4.1
R 


not_normalized_to_save = "NOTnormalized_anina_final_NOV132023.RDS"

normalized_to_save = "normalized_single_and_bulk_final_NOV0923.RDS"

path_to_results_bulk = "q_results_sMRI2/"

path_to_results_single = "/q_results_single/"

cov_files="cov_sMRI_VALUE_AUG102022.RDS" #covariate file

rm(list=ls())
library(data.table)
library(readxl)
library(stringr)


##dont normalize

setwd(path_to_results_bulk)
	# read in data
	vol_bulk <- read.table("aseg.vol.table", header = T)
	lh_thickness_2009_bulk <- read.table("lh.aparc.a2009.thickness.table", header = T)
	rh_thickness_2009_bulk <- read.table("rh.aparc.a2009.thickness.table", header = T)
	lh_area_bulk <- read.table("lh.aparc.area.table", header = T)
	rh_area_bulk <- read.table("rh.aparc.area.table", header = T)
	lh_thickness_bulk <- read.table("lh.aparc.thickness.table", header = T)
	rh_thickness_bulk <- read.table("rh.aparc.thickness.table", header = T)


setwd(path_to_results_single)

	# read in data
	vol_single <- read.table("aseg.vol.table", header = T)
	lh_thickness_2009_single <- read.table("lh.aparc.a2009.thickness.table", header = T)
	rh_thickness_2009_single <- read.table("rh.aparc.a2009.thickness.table", header = T)
	lh_area_single <- read.table("lh.aparc.area.table", header = T)
	rh_area_single <- read.table("rh.aparc.area.table", header = T)
	lh_thickness_single <- read.table("lh.aparc.thickness.table", header = T)
	rh_thickness_single <- read.table("rh.aparc.thickness.table", header = T)

#format
	#vol
		rownames(vol_bulk) <- vol_bulk[,1]
		vol_bulk <- vol_bulk[,-1]
		vol_bulk <- as.data.frame(t(vol_bulk))
		dim(vol_bulk)

		rownames(vol_single) <- vol_single[,1]
		vol_single <- vol_single[,-1]
		vol_single <- as.data.frame(t(vol_single))
		dim(vol_single)

		vol <- cbind(vol_single, vol_bulk)
		
		dim(vol)

	#lh_thickness_2009

		rownames(lh_thickness_2009_bulk) <- lh_thickness_2009_bulk[,1]
		lh_thickness_2009_bulk <- lh_thickness_2009_bulk[,-1]
		lh_thickness_2009_bulk <- as.data.frame(t(lh_thickness_2009_bulk))
		dim(lh_thickness_2009_bulk)

		rownames(lh_thickness_2009_single) <- lh_thickness_2009_single[,1]
		lh_thickness_2009_single <- lh_thickness_2009_single[,-1]
		lh_thickness_2009_single <- as.data.frame(t(lh_thickness_2009_single))
		dim(lh_thickness_2009_single)
		

		lh_thickness_2009 <- cbind(lh_thickness_2009_single, lh_thickness_2009_bulk)
		dim(lh_thickness_2009)

	#rh_thickness_2009
		rownames(rh_thickness_2009_bulk) <- rh_thickness_2009_bulk[,1]
		rh_thickness_2009_bulk <- rh_thickness_2009_bulk[,-1]
		rh_thickness_2009_bulk <- as.data.frame(t(rh_thickness_2009_bulk))
		dim(rh_thickness_2009_bulk)

		rownames(rh_thickness_2009_single) <- rh_thickness_2009_single[,1]
		rh_thickness_2009_single <- rh_thickness_2009_single[,-1]
		rh_thickness_2009_single <- as.data.frame(t(rh_thickness_2009_single))
		dim(rh_thickness_2009_single)

		rh_thickness_2009 <- cbind(rh_thickness_2009_single, rh_thickness_2009_bulk)
		dim(rh_thickness_2009)

	#lh_area
		rownames(lh_area_bulk) <- lh_area_bulk[,1]
		lh_area_bulk <- lh_area_bulk[,-1]
		lh_area_bulk <- as.data.frame(t(lh_area_bulk))
		dim(lh_area_bulk)

		rownames(lh_area_single) <- lh_area_single[,1]
		lh_area_single <- lh_area_single[,-1]
		lh_area_single <- as.data.frame(t(lh_area_single))
		dim(lh_area_single)

		lh_area <- cbind(lh_area_single, lh_area_bulk)
		dim(lh_area)

	#rh_area
		rownames(rh_area_bulk) <- rh_area_bulk[,1]
		rh_area_bulk <- rh_area_bulk[,-1]
		rh_area_bulk <- as.data.frame(t(rh_area_bulk))
		dim(rh_area_bulk)

		rownames(rh_area_single) <- rh_area_single[,1]
		rh_area_single <- rh_area_single[,-1]
		rh_area_single <- as.data.frame(t(rh_area_single))
		dim(rh_area_single)

		rh_area <- cbind(rh_area_single, rh_area_bulk)
		dim(rh_area)

	#lh_thickness
		rownames(lh_thickness_bulk) <- lh_thickness_bulk[,1]
		lh_thickness_bulk <- lh_thickness_bulk[,-1]
		lh_thickness_bulk <- as.data.frame(t(lh_thickness_bulk))
		dim(lh_thickness_bulk)

		rownames(lh_thickness_single) <- lh_thickness_single[,1]
		lh_thickness_single <- lh_thickness_single[,-1]
		lh_thickness_single <- as.data.frame(t(lh_thickness_single))
		dim(lh_thickness_single)

		lh_thickness <- cbind(lh_thickness_single, lh_thickness_bulk)
		dim(lh_thickness)

	#rh_thickness
		rownames(rh_thickness_bulk) <- rh_thickness_bulk[,1]
		rh_thickness_bulk <- rh_thickness_bulk[,-1]
		rh_thickness_bulk <- as.data.frame(t(rh_thickness_bulk))
		dim(rh_thickness_bulk)

		rownames(rh_thickness_single) <- rh_thickness_single[,1]
		rh_thickness_single <- rh_thickness_single[,-1]
		rh_thickness_single <- as.data.frame(t(rh_thickness_single))
		dim(rh_thickness_single)

		rh_thickness <- cbind(rh_thickness_single, rh_thickness_bulk)
		dim(rh_thickness)


	rownames(lh_area) <- gsub("\\.", "_", rownames(lh_area))
	rownames(rh_thickness) <- gsub("\\.", "_", rownames(rh_thickness))
	rownames(rh_area) <- gsub("\\.", "_", rownames(rh_area))

	rownames(vol) <- gsub("X4th.Ventricle", "Fourth_Ventricle", rownames(vol))
	rownames(vol) <- gsub("X3rd.Ventricle", "Third_Ventricle", rownames(vol))
	rownames(vol) <- gsub("X5th.Ventricle", "Fifth_Ventricle", rownames(vol))
	rownames(vol) <- gsub("Right_Accumbens_area", "right_accumbens", rownames(vol))
	rownames(vol) <- gsub("Left_Accumbens_area", "left_accumbens", rownames(vol))
	rownames(vol) <- gsub("\\.", "_", rownames(vol))

	norm1 <- rbind(vol,lh_area,lh_thickness, rh_area,rh_thickness) 
	norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent1",]
	norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent2",]
	norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent3",]
	norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent4",]
	norm1=norm1[!rownames(norm1)=="eTIV1",]
	norm1=norm1[!rownames(norm1)=="eTIV2",]
	norm1=norm1[!rownames(norm1)=="eTIV3",]	
	

##unnormalized data:

setwd(path_to_results_bulk)
	# read in data
	vol_bulk <- read.table("aseg.vol.table", header = T)
	lh_thickness_2009_bulk <- read.table("lh.aparc.a2009.thickness.table", header = T)
	rh_thickness_2009_bulk <- read.table("rh.aparc.a2009.thickness.table", header = T)
	lh_area_bulk <- read.table("lh.aparc.area.table", header = T)
	rh_area_bulk <- read.table("rh.aparc.area.table", header = T)
	lh_thickness_bulk <- read.table("lh.aparc.thickness.table", header = T)
	rh_thickness_bulk <- read.table("rh.aparc.thickness.table", header = T)


setwd(path_to_results_single)

	# read in data
	vol_single <- read.table("aseg.vol.table", header = T)
	lh_thickness_2009_single <- read.table("lh.aparc.a2009.thickness.table", header = T)
	rh_thickness_2009_single <- read.table("rh.aparc.a2009.thickness.table", header = T)
	lh_area_single <- read.table("lh.aparc.area.table", header = T)
	rh_area_single <- read.table("rh.aparc.area.table", header = T)
	lh_thickness_single <- read.table("lh.aparc.thickness.table", header = T)
	rh_thickness_single <- read.table("rh.aparc.thickness.table", header = T)

##format
	
	#vol
		rownames(vol_bulk) <- vol_bulk[,1]
		vol_bulk <- vol_bulk[,-1]
		vol_bulk <- as.data.frame(t(vol_bulk))
		dim(vol_bulk)

		rownames(vol_single) <- vol_single[,1]
		vol_single <- vol_single[,-1]
		vol_single <- as.data.frame(t(vol_single))
		dim(vol_single)

		vol <- cbind(vol_single, vol_bulk)
		
		dim(vol)

	#lh_thickness_2009
		rownames(lh_thickness_2009_bulk) <- lh_thickness_2009_bulk[,1]
		lh_thickness_2009_bulk <- lh_thickness_2009_bulk[,-1]
		lh_thickness_2009_bulk <- as.data.frame(t(lh_thickness_2009_bulk))
		dim(lh_thickness_2009_bulk)

		rownames(lh_thickness_2009_single) <- lh_thickness_2009_single[,1]
		lh_thickness_2009_single <- lh_thickness_2009_single[,-1]
		lh_thickness_2009_single <- as.data.frame(t(lh_thickness_2009_single))
		dim(lh_thickness_2009_single)

		lh_thickness_2009 <- cbind(lh_thickness_2009_single, lh_thickness_2009_bulk)
		dim(lh_thickness_2009)

	#rh_thickness_2009
		rownames(rh_thickness_2009_bulk) <- rh_thickness_2009_bulk[,1]
		rh_thickness_2009_bulk <- rh_thickness_2009_bulk[,-1]
		rh_thickness_2009_bulk <- as.data.frame(t(rh_thickness_2009_bulk))
		dim(rh_thickness_2009_bulk)

		rownames(rh_thickness_2009_single) <- rh_thickness_2009_single[,1]
		rh_thickness_2009_single <- rh_thickness_2009_single[,-1]
		rh_thickness_2009_single <- as.data.frame(t(rh_thickness_2009_single))
		dim(rh_thickness_2009_single)

		rh_thickness_2009 <- cbind(rh_thickness_2009_single, rh_thickness_2009_bulk)
		dim(rh_thickness_2009)


	#lh_area
		rownames(lh_area_bulk) <- lh_area_bulk[,1]
		lh_area_bulk <- lh_area_bulk[,-1]
		lh_area_bulk <- as.data.frame(t(lh_area_bulk))
		dim(lh_area_bulk)

		rownames(lh_area_single) <- lh_area_single[,1]
		lh_area_single <- lh_area_single[,-1]
		lh_area_single <- as.data.frame(t(lh_area_single))
		dim(lh_area_single)

		lh_area <- cbind(lh_area_single, lh_area_bulk)
		dim(lh_area)

	#rh_area
		rownames(rh_area_bulk) <- rh_area_bulk[,1]
		rh_area_bulk <- rh_area_bulk[,-1]
		rh_area_bulk <- as.data.frame(t(rh_area_bulk))
		dim(rh_area_bulk)

		rownames(rh_area_single) <- rh_area_single[,1]
		rh_area_single <- rh_area_single[,-1]
		rh_area_single <- as.data.frame(t(rh_area_single))
		dim(rh_area_single)

		rh_area <- cbind(rh_area_single, rh_area_bulk)
		dim(rh_area)

	#lh_thickness
		rownames(lh_thickness_bulk) <- lh_thickness_bulk[,1]
		lh_thickness_bulk <- lh_thickness_bulk[,-1]
		lh_thickness_bulk <- as.data.frame(t(lh_thickness_bulk))
		dim(lh_thickness_bulk)

		rownames(lh_thickness_single) <- lh_thickness_single[,1]
		lh_thickness_single <- lh_thickness_single[,-1]
		lh_thickness_single <- as.data.frame(t(lh_thickness_single))
		dim(lh_thickness_single)

		lh_thickness <- cbind(lh_thickness_single, lh_thickness_bulk)
		dim(lh_thickness)

	#rh_thickness
		rownames(rh_thickness_bulk) <- rh_thickness_bulk[,1]
		rh_thickness_bulk <- rh_thickness_bulk[,-1]
		rh_thickness_bulk <- as.data.frame(t(rh_thickness_bulk))
		dim(rh_thickness_bulk)

		rownames(rh_thickness_single) <- rh_thickness_single[,1]
		rh_thickness_single <- rh_thickness_single[,-1]
		rh_thickness_single <- as.data.frame(t(rh_thickness_single))
		dim(rh_thickness_single)

		rh_thickness <- cbind(rh_thickness_single, rh_thickness_bulk)
		dim(rh_thickness)


##lh_area i.e. lh.aparc.area.table:
				#value/lh_WhiteSurfArea_area
				for(i in 1:ncol(lh_area)) {       
				  lh_area[1:34,i] <- lh_area[1:34,i] / lh_area[35,i]
				}

		rownames(lh_area) <- gsub("\\.", "_", rownames(lh_area))

		#lh_thickness i.e lh.aparc.thickness 
				#value/lh_MeanThickness_thickness
				for(i in 1:ncol(lh_thickness)) {       
				  lh_thickness[1:34,i] <- lh_thickness[1:34,i] / lh_thickness[35,i]
				}

# rh_thickness i.e rh.aparc.thickness
				#value/rh_MeanThickness_thickness
				for(i in 1:ncol(rh_thickness)) {       
				  rh_thickness[1:34,i] <- rh_thickness[1:34,i] / rh_thickness[35,i]
				}

		rownames(rh_thickness) <- gsub("\\.", "_", rownames(rh_thickness))


		#rh_area i.e rh.aparc.area.table
				#value/rh_WhiteSurfArea_area
				for(i in 1:ncol(rh_area)) {       
				  rh_area[1:34,i] <- rh_area[1:34,i] / rh_area[35,i]
				}

		rownames(rh_area) <- gsub("\\.", "_", rownames(rh_area))

		for(i in 1:ncol(vol)) {       
				  vol[1:60,i] <- vol[1:60,i] / vol[66,i]
				}
		rownames(vol) <- gsub("X4th.Ventricle", "Fourth_Ventricle", rownames(vol))
		rownames(vol) <- gsub("X3rd.Ventricle", "Third_Ventricle", rownames(vol))
		rownames(vol) <- gsub("X5th.Ventricle", "Fifth_Ventricle", rownames(vol))
	    rownames(vol) <- gsub("Right_Accumbens_area", "right_accumbens", rownames(vol))
	    rownames(vol) <- gsub("Left_Accumbens_area", "left_accumbens", rownames(vol))
		rownames(vol) <- gsub("\\.", "_", rownames(vol))

		norm1 <- rbind(vol,lh_area,lh_thickness, rh_area,rh_thickness) 

		norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent1",]
		norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent2",]
		norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent3",]
		norm1=norm1[!rownames(norm1)=="BrainSegVolNotVent4",]
		norm1=norm1[!rownames(norm1)=="eTIV1",]
		norm1=norm1[!rownames(norm1)=="eTIV2",]
		norm1=norm1[!rownames(norm1)=="eTIV3",]

		dim(norm1)

saveRDS(norm1, normalized_to_save)


