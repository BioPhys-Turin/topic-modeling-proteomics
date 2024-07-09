from nsbm import nsbm
import pandas as pd
import numpy as np
import pickle 
import sys, os

#this is for the single omics analysis using the nSBM implementation

#medium
rna = pd.read_csv("data/sbm_input_brca_cptac/brca_rna_119s_hsbmprep_4k.csv", index_col=0).fillna(0).astype(int)
prot = pd.read_csv("data/sbm_input_brca_cptac/brca_prot_119s_hsbmprep_2k.csv", index_col=0).fillna(0).astype(int)
phos = pd.read_csv("data/sbm_input_brca_cptac/brca_phos_119s_hsbmprep_8k.csv", index_col=0).fillna(0).astype(int)


nsbm_model = nsbm()
nsbm_model.make_graph_multiple_df(rna, [prot, phos])
nsbm_model.fit()

os.chdir("sbm_output_brca_cptac")
file = open('brca_nsbm_medium', 'wb')
pickle.dump(nsbm_model, file)
file.close()

os.chdir("brca_nsbm_tests_medium")

nsbm_model.save_data()

#here saving the worddist file is not necessary as it doen it automatically

#you can repeat this block for all the different omics layers adapted to your filepaths