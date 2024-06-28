# Using stochastic based modelling (SBM) topic modelling on proteomics and phosphoproteomics data (single-omics & multi-omics integration)

In this is a tutorial and application example for the hSBM algorithm (single-omics) by Martin Gerlach 
https://github.com/martingerlach/hSBM_Topicmodel 
and the nSBM implementation (multi-omics) by the Biophys Turin group 
[https://github.com/BioPhys-Turin/nsbm.](https://github.com/martingerlach/hSBM_Topicmodel )


Gerlach, M., Peixoto, T. P., & Altmann, E. G. (2018). A network approach to topic models. Science Advances, 4(7), eaaq1360. https://doi.org/10.1126/sciadv.aaq1360

Get second citation

Data used in this work is the pulic CPTAC data set available via the python package cptac https://pypi.org/project/cptac/, it includes bulk transcriptomics,
proteomics & phosphoproteomics datasets as well as extensice clinical data for 11 cancer types. In the code here only uses the Breast Cancer cohort 
(119 patients), but it can be applied to all other cohorts. Further clinical information for each cancer type cohort has been published can be used as further ressource, 
like in this case for the breast cancer (https://www.cell.com/cell/fulltext/S0092-8674(20)31400-8).

## Structure of the tutorial

1. Preprocessing (R markdown)
2. Run SBM analysis, nSBM & hSBM (Python, needs to be done on computing cluster, min of 30gb ram required)
3. Readout of all results, visualisation (R markdown)
4. Further investigation of results (R markdown)

The results are provided, so if you don't want/ can't run the analysis, you can directly start with the script 3 without having to install anything additional to the R packages specified in the markdown.
Still, make sure to check out script 1 with preprocessing script first, to properly understand the input data and matching results.

## Set Up SBM

If you want to run SBM on your own data, first follow the setup requirements provided by the hSBM tutorial https://github.com/martingerlach/hSBM_Topicmodel. 
Then install the nsbm python package from the nSBM tutorial https://github.com/martingerlach/hSBM_Topicmodel. 

```
conda create --name sbm_env python=3.7
conda activate sbm_env
conda install -c conda-forge gtk3 pygobject matplotlib graph-tool #hSBM requirements
conda install -c conda-forge nsbm #nsbm requirements

```

Make sure to additionally download the hSBM code base from here `git clone https://github.com/martingerlach/hSBM_Topicmodel.git` 
(theoretically it is contained in the nsbm python package but i haven't tested it so lets just be safe). 

R requirements ? should I add a txt with all the info?







```

