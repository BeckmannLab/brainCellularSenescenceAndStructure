
#scenic run 
cell=(All) #run also for Exc and MG
cell=(Exc)
cell=(MG)

DIR=/SCENIC_working

CELLDIR=${DIR}/ssAnina
exprMat=${CELLDIR}/${cell}_exprDat_Final_ssAnina.csv
adjOut=${CELLDIR}/${cell}_adj_Final_ssAnina.csv
regOut=${CELLDIR}/${cell}_reg_Final_ssAnina.csv
aucOut=${CELLDIR}/${cell}_auc_Final_ssAnina.csv

module load pyscenic

pyscenicVenv

arboreto_with_multiprocessing.py \
    ${exprMat} \
    ${DIR}/hs_hgnc_tfs.txt \
    --method grnboost2 \
    --output ${adjOut} \
    --num_workers 30 \
    --seed 777 \
    --transpose

pyscenic ctx ${adjOut} \
    ${DIR}/*feather \
    --annotations_fname ${DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tsv \
    --expression_mtx_fname ${exprMat} \
    --output ${regOut} \
    --mask_dropouts \
    --num_workers 30 \
    --transpose

pyscenic aucell \
    ${exprMat} \
    ${regOut} \
    --output ${aucOut} \
    --num_workers 30 \
    --transpose
