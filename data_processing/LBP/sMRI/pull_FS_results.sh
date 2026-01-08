#pull quant from fs results in bash 	
	cp -r output_fs output_fs2

	outputDIR=output_fs2
	output_path=q_results_sMRI2


	#run the following (order matters and do it one at a time. copying and pasting whole thing seems    to error)
	module load freesurfer
	DIR=outputDIR
	mkdir q_results_sMRI2 #output folder
	cd ${DIR}
	SUBJECTS_DIR=outputDIR
	for d in `ls -d *`; do
	    aparcstats2table --hemi lh --subjects scan* --meas thickness --parc aparc.a2009s --tablefile output_path/lh.aparc.a2009.thickness.table
		#This is right hemisphere thickness table
		aparcstats2table --hemi rh --subjects scan* --meas thickness --parc aparc.a2009s --tablefile output_path/rh.aparc.a2009.thickness.table
		#This is left hemisphere
		aparcstats2table --hemi lh --subjects scan* --meas thickness --tablefile output_path/lh.aparc.thickness.table
		#This is right hemisphere
		aparcstats2table --hemi rh --subjects scan* --meas thickness --tablefile output_path/rh.aparc.thickness.table
		#Surface Area Left hemisphere
		aparcstats2table --hemi lh --subjects scan* --tablefile output_path/lh.aparc.area.table
		#Surface Area right hemisphere
		aparcstats2table --hemi rh --subjects scan* --tablefile output_path/rh.aparc.area.table
		#Subcortical Volume
		asegstats2table --subjects scan* --tablefile output_path/aseg.vol.table
		#Mean intensity
		asegstats2table --subjects scan* --meas mean --tablefile output_path/aseg.mean-intensity.table    
	done 
	
