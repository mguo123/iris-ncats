### MARGARET GUO
### CREATED ON: 11/02/17 
### IMPORT ALL MODULES AND FILES NEEDED FOR the NCATS module

import sys, os
import csv
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import pickle
import networkx as nx
import numpy as np
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from textwrap import wrap
from collections import defaultdict
import itertools, time

# to import the rest of the modules more easily
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) # points to scripts DIR

import run_main
from get_results import get_results
from get_output_data import get_output_data
import find_neighborhood_beta
from find_neighborhood_beta import find_neighborhood as fgn
import run_analysis as ann
from get_associations_deprecated import get_associations

# import association file data
# Define variables with paths
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # should yield the same as .. (with subfolders rcsc, scripts, and results)
RSCS_DIR = os.path.join(PARENT_DIR, 'rscs')
RESULTS_DIR = os.path.join(PARENT_DIR, 'results')


# loading data: interactome
NETWORKF = os.path.join(RSCS_DIR,'merged_interact_netx.pkl')#non-specific interaction network 

print('loading interactome data from', NETWORKF)
GENE_GRAPH = pickle.load(open(NETWORKF,'rb'))
print("(((((((((()))))))))")
print(GENE_GRAPH.nodes())
raise ValueError("Passed!!")

# loading random data: 
ann.rand_dir = os.path.join(RESULTS_DIR, 'rand_iRefplus_intome/summary/')

# load drug targets, then format
DTF = os.path.join(RSCS_DIR, 'drug_intome_targets.pkl')
print('loading drug targets', DTF)
DTD = pickle.load(open(DTF,'rb'))

print('loading graph files')
G_TO_RSIDS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rsid_eQTL.pkl'),'rb'))
RS_TO_DATA = pickle.load(open(os.path.join(RSCS_DIR,'rsid_g_pval_rsqr.pkl'),'rb')) # rsid_g_pval_rsqr[rs]=[g,pv,rsq] 
G_TO_DISGENNET = pickle.load(open(os.path.join(RSCS_DIR,'disGeNet_gene_dis_score_dict.pkl'),'rb')) #
G_TO_OMIM = pickle.load(open(os.path.join(RSCS_DIR,'OMIM_genes_to_phenotypes.pkl'),'rb')) # OMIM phenotypes
G_TO_PHW_SNPS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rs_phWAS.pkl'),'rb')) # gene to SNPs from PheWAS
PHW_SNPS_TO_PHEN = pickle.load(open(os.path.join(RSCS_DIR,'rs_to_phenOdds_phWAS.pkl'),'rb')) # SNPs to phenotype
ALL_ASSOC = pickle.load(open(os.path.join(RSCS_DIR,'all_assoc_to_nodes.pkl'),'rb'))# all associations, and their genes/SNPs
INTOME_SZIE = pickle.load(open(os.path.join(RSCS_DIR,'interactome_size.pkl'),'rb'))

# Analysis parameters

SCORE_THRESHOLD_START = 0.8
SCORE_THRESHOLD_MIN = 0.5
SCORE_DELTA = 0.05