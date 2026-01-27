# brainCellularSenesenceAndStructure
This repository contains code used in Lund et al., *Establishing the relationship between brain cellular senescence and brain structure*. This repository is organized to reflect analyses reported in the manuscript

### Data Availability

This study includes both newly generated data and data obtained from previously published or publicly available datasets.

#### LBP

All primary datasets generated as part of this study, along with associated metadata, are publicly available through Synapse.org under the study Synapse project **syn26337520**. These include:

- **Single-nucleus RNA-seq (snRNA-seq)**: syn65929068, syn65929069  
- **Bulk RNA-seq**: syn27127036  
- **Mass spectrometry (MS) proteomics**: syn65929055  
- **Structural MRI (sMRI) features**:  
  - Raw features: syn65929054  
  - Summarized features: syn65929053
  - Associated metadata: syn27127033

For access to raw images to replicate study results, please reach out directly to corresponding authors. Raw images are stored on the Minerva computing cluster at the Icahn School of Medicine at Mount Sinai; access requires institutional IRB approval.

#### External Datasets Used in This Study

This study also incorporates data from previously published cohorts and public resources, subject to their respective access requirements:

- **brainSCOPE**  
  - Counts: https://brainscope.gersteinlab.org/output-sample-annotated-matrix.html  
  - Metadata:  
    - https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt  
    - https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_mapping.xlsx  

- **Neurodevelopmental dataset (Herring et al.)**  
  - GEO: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168408/suppl/filelist.txt  
  - Associated publication: https://www.cell.com/cell/pdf/S0092-8674(22)01258-2.pdf  

- **Tang et al.**: GEO GSE115301  
- **Chan et al.**: GEO GSE175533  
- **Teo et al.**: GEO GSE119807  

- **GTEx v8** (accessed June 5, 2017):  
  https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression  

- **ROSMAP**: Synapse syn4164376  
- **MSBB**: Synapse syn7416949  


### Reproducibility Notes

- Scripts in `senescence_pipeline/` provide a general example of the senescence definition workflow.
- Dataset-specific scripts reproduce analyses as performed in the manuscript but may require dataset-specific preprocessing, metadata, and access permissions.
- Paths, compute environments, and runtime parameters may need adjustment for reuse.

```text
brainCellularSenescenceAndStructure/
├─ senescence_pipeline/                         general example script for defining senescence
├─ data_processing/                         Preprocessing of data
│  ├─ BrainSCOPE/                           Data preprocessing, defining senescence, pseudobulk
│  ├─ Herring/                              Data preprocessing, defining senescence, pseudobulk
│  ├─ LBP/sMRI/                             sMRI data preprocessing
│  ├─ MSBB/                                 MSBB data alignment to genome reference (RAPiD)
│  └─ ROSMAP/                               ROSMAP data alignment to genome reference (RAPiD)
│
├─ QC/                                      QC analyses
│  ├─ senCID/                               senCID analyses in bulk and snRNA-seq
│  ├─ gtex_gene_expr_comparison.R           Code for GTEx brain region gene expression
│  ├─ multi-tissue_PFC_comparison.R         Code for PFC gene expression across brain and tissues
│  └─ sen_proportion_vs_sMRI.R              Code for senescent cell proportion vs sMRI features
│
├─ DE_analyses/                             Differential expression analyses
│  ├─ Binary_DE/                            Binary differential expression analyses
│  │  ├─ BrainSCOPE/                        DE analysis for BrainSCOPE data
│  │  ├─ Chan/                              defining senescence, pseudobulk, and DE analysis for Chan et al.
│  │  ├─ Herring/                           DE analysis for Herring et al. data
│  │  ├─ LBP/                               defining senescence, pseudobulk, and DE analysis for LBP data (biopsy and bank)
│  │  ├─ Tang/                              defining senescence, pseudobulk, and DE analysis for Tang et al.
│  │  └─ Teo/                               defining senescence, pseudobulk, and DE analysis for Teo et al.
│  │
│  └─ Continuous_DE/                        Continuous differential expression analyses
│     ├─ MS_sMRI_DE/                        DE analysis for mass spectrometry and sMRI
│     ├─ bulkRNAseq_sMRI_DE/                DE analysis for bulk RNA-seq and sMRI
│     └─ snRNAseq_sMRI_DE/                  DE analysis for snRNA-seq and sMRI
│
├─ Pathway_enrichment/                      Pathway enrichment analyses from DE results
│
├─ SCENIC/                                  Infer gene regulatory networks
│
├─ Figures/                                 Figure generation scripts
│  ├─ Main/                                 Scripts for main manuscript figures
│  └─ Supplementary/                        Scripts for supplementary figures
│
├─ functions/                               Functions used in senescence analyses
│
└─ misc/                                    mapping and phenotype files
```

### Contact

Questions regarding data access or analyses should be directed to anina.lund@mssm.edu

### Citation

If you use this code, please cite:

Lund, A. N., et al. (2026). *Establishing the relationship between brain cellular senescence and brain structure*. **Cell**. https://doi.org/10.1016/j.cell.2025.10.014
