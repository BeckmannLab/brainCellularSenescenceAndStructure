`defining_senescence_pipeline_example.R` is a general example script for defining senescence in single-cell and single-nuclei RNA-seq data. It illustrates the core analysis workflow used across datasets.

This script requires the following helper functions:

- `functions/dreamletCompareClisters_edgeR.R`
- `functions/processOneAssay_edgeR.R`

The following scripts were used in the paper to define senescence in specific datasets. These scripts adapt the general pipeline to individual datasets and perform dataset-specific senescence definition and differential expression (DE) analyses:
- `DE_analyses/Binary_DE/BrainSCOPE/DE_jobs_bs.R` BrainSCOPE
- `DE_analyses/Binary_DE/Chan/Chan_validation.R` Chan et al
- `DE_analyses/Binary_DE/Herring/Herring_DE.R` Herring et al (neurodevelopment)
- `DE_analyses/Binary_DE/LBP/liv_sen_DE.R` LBP biopsy
- `DE_analyses/Binary_DE/LBP/pm_sen_DE.R` LBP bank
- `DE_analyses/Binary_DE/Tang/Tang_validation_DE.R` Tang et al
- `DE_analyses/Binary_DE/Teo/Teo_validation.R` Teo et al
