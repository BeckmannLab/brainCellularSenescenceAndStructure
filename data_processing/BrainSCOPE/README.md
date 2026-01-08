## Preprocessing BrainSCOPE

To preprocess **BrainSCOPE**, run the scripts in the following order:

1. **Set up individual datasets**  
   Folder: `set_up/`  

2. **Subset by cell type**  
   Folder: `subset_celltype/`  

3. **Run AUCell by cell type**  
   Script: `aucell/aucell_by_celltype.sh`  

4. **Define thresholds across datasets**  
   Script: `aucell/across_ds_thresholds.R`  

5. **Apply thresholds per dataset**  
   Folder: `aucell/`  

6. **Convert to pseudobulk**  
   Folder: `pseudobulk/`  

7. **Combine pseudobulk datasets**  
   Script: `pseudobulk/combine_pb.R`  
