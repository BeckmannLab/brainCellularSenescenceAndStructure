##this made the scan paths and folders so we can run freesurfer 
DIR=... #path to where the DICOM files are
cd ${DIR}
mkdir ${DIR}/output 
while read LINE
do 
CHECK=`echo ${LINE} | awk -F"|" '{print $1}'`
if [[ ${CHECK} != '"path"' ]]
then 
  NEWPATH=${DIR}/output/`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
  OLDPATH=`echo ${LINE} | awk -F "|" '{print $1}'`
  if [[ ! -d ${NEWPATH} ]]
  then mkdir ${NEWPATH}
  fi
  eval ln -f -s ${OLDPATH} ${DIR}/input/`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
fi 
done < map_input_output_imaging_sMRI.txt #file with paths to image files

##convert DICOM to nifti
while read LINE
do 
CHECK=`echo ${LINE} | awk -F"|" '{print $1}'`
if [[ ${CHECK} != '"path"' ]]
then 
  SCAN=`echo ${LINE} | awk -F "|" '{print $2}' | sed s/'"'/''/g`
  I=${DIR}/input/${SCAN}
  O=${DIR}/output/${SCAN}
  dcm2niix -o ${O} -f "%f" -p y -z y ${I}
  STATUS=$?   
  if [[ ${STATUS} -eq 0 ]]
  then touch ${O}/${SCAN}.success
  else touch ${O}/${SCAN}.fail;exit 1
  fi
fi 
done < map_input_output_imaging_sMRI.txt #file with paths to image files

#checks
ls ${DIR}/output/scan*/* | grep ".success" | wc -l
ls ${DIR}/output/scan*/* | grep ".fail" | wc -l
ls | wc -l
ls ../input/ | wc -l
