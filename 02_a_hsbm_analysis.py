import pandas as pd
import sys, os
sys.path.append("hSBM_Topicmodel/")
from sbmtm import sbmtm
import pickle

#this is for the single omics analysis using the hSBM code

rna = pd.read_csv("data/sbm_input_brca_cptac/brca_rna_119s_hsbmprep_4k.csv", index_col=0).fillna(0).astype(int)
prot = pd.read_csv("data/sbm_input_brca_cptac/brca_prot_119s_hsbmprep_2k.csv", index_col=0).fillna(0).astype(int)
phos = pd.read_csv("data/sbm_input_brca_cptac/brca_phos_119s_hsbmprep_8k.csv", index_col=0).fillna(0).astype(int)


#define model
model = sbmtm()
#Build graph from one omics layer
model.make_graph_from_BoW_df(df=rna)
#fit the model (has to be done on a computing cluster, only protc)
model.fit()

os.chdir("sbm_output_brca_cptac") # the directory for all your sbm results
file = open('brca_nsbm_medium', 'wb') #the file to save the entire model
pickle.dump(model, file)
file.close()

os.chdir("brca_rna_tests_4k") #file to save all your matrices, create the file beforehand or replace this line with os.mkdir() command
model.save_data() #save all the matrices in the file

for i in range(0,4): #extract the worddist matrix manually as it doensn't get extracted with save_data()
    results = model.get_groups(l=i)
    df = pandas.DataFrame(data=results['p_w_tw'])
    df = df.assign(word_names=model.words)
    f = 'brca_prot_4k_l' + str(i) +'_worddist.csv'
    df.to_csv(f, index=False)


os.chdir("..")

#you can repeat this block for all the different omics layers adapted to your filepaths