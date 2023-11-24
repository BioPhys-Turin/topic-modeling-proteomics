from nsbm import nsbm
import pandas as pd
import numpy as np
import sys, os
import pickle



brca_rna = pd.read_csv("nsbm_input/brca_rna_119s_nsbmprep.csv", index_col=0).fillna(0).astype(int)
brca_rna = brca_rna/ 100
brca_prot = pd.read_csv("nsbm_input/brca_proteo_119s_nsbmprep.csv", index_col=0).fillna(0).astype(int)
brca_prot = brca_prot/ 1000
brca_phos = pd.read_csv("nsbm_input/brca_phospho_119s_nsbmprep.csv", index_col=0).fillna(0).astype(int)
brca_phos = brca_phos/ 1000

nsbm_model = nsbm()
nsbm_model.make_graph_multiple_df(brca_rna, [brca_prot, brca_phos])
nsbm_model.fit()

#os.chdir("nsbm_output")
#file = open('brca_nsbm_3layers', 'wb')
#pickle.dump(nsbm_model, file)
#file.close()

#os.chdir("brca_3layers")

#nsbm_model.save_data()
