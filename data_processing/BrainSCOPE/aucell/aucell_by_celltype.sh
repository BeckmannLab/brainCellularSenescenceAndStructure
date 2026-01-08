#in bash combine files for each cell type across ds

##by cell type
cat Exc*.csv  > all_Exc.csv

cat Micro*.csv  > all_Micro.csv

cat OPC*.csv > all_OPC.csv

cat Astro*.csv > all_Astro.csv

cat Int*.csv > all_Int.csv

cat Oligo*.csv > all_Oligo.csv

#############
#python
#############

#set up
module load anaconda3/2021.5
conda activate pyscenic

ml python/3.10
python

import pandas as pd
import os
import glob
import scanpy as sc
import scipy.sparse as sp
import numpy as np
import pandas as pd
import ctxcore
from ctxcore.genesig import GeneSignature
from pyscenic.aucell import create_rankings, enrichment
from pyscenic.aucell import aucell, derive_auc_threshold, create_rankings
import anndata as ad
import os
import re
import polars as pl
from scipy.sparse import csr_matrix

#######################
#get expression object
#######################

####################
#set up geneset file
##################### 
gs= pd.read_csv('signature_146.tsv', sep = "\t") #gene list 

gene_list = gs['signature'].dropna().tolist()  # Remove any NaN values if they exist
gene_set_name = "MyGeneSet"  # You can choose any name for your gene set
gene_set = {gene_set_name: set(gene_list)}

# Convert the gene set to a GeneSignature object
gene_signatures = [GeneSignature(name=name, gene2weight=dict(zip(genes, [1]*len(genes)))) 
                   for name, genes in gene_set.items()]

#############
#aucell
############

######exc

df_exc = pl.read_csv("all_Exc.csv", ignore_errors=True,truncate_ragged_lines=True)

df_panda = pd.DataFrame(df_exc[:,1:].to_numpy().astype('float32'))

df_panda.index = df_exc[:, 0].to_list()

# Assign column names (excluding the first column)
df_panda.columns = df_exc.columns[1:]

# df_panda.max() - df_panda.min()
del df_exc  
aucs_mtx_exc = aucell(df_panda, gene_signatures, num_workers=50)
aucs_mtx_exc.to_csv('exc_aucell_across.csv', index=True) 

#############
#micro
############
df_micro = pl.read_csv("all_Micro.csv", ignore_errors=True,truncate_ragged_lines=True)

df_panda = pd.DataFrame(df_micro[:,1:].to_numpy().astype('float32'))


names = df_micro[:, 0].to_list()
"Micro_1041" in names

df_panda.index = df_micro[:, 0].to_list()

# Assign column names (excluding the first column)
df_panda.columns = df_micro.columns[1:]


# df_panda.max() - df_panda.min()
del df_micro
aucs_mtx_micro = aucell(df_panda, gene_signatures, num_workers=50)
aucs_mtx_micro.to_csv('micro_aucell_across.csv', index=True)


#############
#Int
############
df_int = pl.read_csv("all_Int.csv", ignore_errors=True,truncate_ragged_lines=True)

df_panda = pd.DataFrame(df_int[:,1:].to_numpy().astype('float32'))

df_panda.index = df_int[:, 0].to_list()

# Assign column names (intluding the first column)
df_panda.columns = df_int.columns[1:]

# df_panda.max() - df_panda.min()
del df_int  
aucs_mtx_int = aucell(df_panda, gene_signatures, num_workers=50)
aucs_mtx_int.to_csv('int_aucell_across.csv', index=True) 


#############
#oli
############
df_oli = pl.read_csv("all_Oligo.csv", ignore_errors=True,truncate_ragged_lines=True)
df_panda = pd.DataFrame(df_oli[:,1:].to_numpy().astype('float32'))

df_panda.index = df_oli[:, 0].to_list()

# Assign column names (intluding the first column)
df_panda.columns = df_oli.columns[1:]

# df_panda.max() - df_panda.min()
del df_oli  
aucs_mtx_oli = aucell(df_panda, gene_signatures, num_workers=50)

aucs_mtx_oli.to_csv('oli_aucell_across.csv', index=True) 

#############
#OPC
############
df_OPC = pl.read_csv("all_OPC.csv", ignore_errors=True,truncate_ragged_lines=True)

df_panda = pd.DataFrame(df_OPC[:,1:].to_numpy().astype('float32'))

df_panda.index = df_OPC[:, 0].to_list()

# Assign column names (intluding the first column)
df_panda.columns = df_OPC.columns[1:]

# df_panda.max() - df_panda.min()
del df_OPC  
aucs_mtx_OPC = aucell(df_panda, gene_signatures, num_workers=50)

aucs_mtx_OPC.to_csv('OPC_aucell_across.csv', index=True) 

#############
#astro
############
df_astro = pl.read_csv("all_Astro.csv", ignore_errors=True,truncate_ragged_lines=True)

df_panda = pd.DataFrame(df_astro[:,1:].to_numpy().astype('float32'))

df_panda.index = df_astro[:, 0].to_list()

# Assign column names (intluding the first column)
df_panda.columns = df_astro.columns[1:]

# df_panda.max() - df_panda.min()
del df_astro  
aucs_mtx_astro = aucell(df_panda, gene_signatures, num_workers=50)

aucs_mtx_astro.to_csv('astro_aucell_across.csv', index=True) 

