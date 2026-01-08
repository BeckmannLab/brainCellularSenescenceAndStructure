
DIR=... #path to where the files are
##run freesurfer
module load freesurfer/6.0.0
mkdir ${DIR}/output_fs/
cd ${DIR}/output/
for d in `ls -d *`; do
        cd ${d}
        bsub -J ${d} -n 4 -P acc_psychgen -W 100:00 -q premium -R "rusage[mem=10000]" -oo "fs.stdout" -eo "fs.stderr" recon-all -i ${d}.nii.gz -s ${d} -sd ${DIR}/output_fs/ -all -qcache;
        cd ..
done

grep Success ${DIR}/output/*.sdtout | wc -l

