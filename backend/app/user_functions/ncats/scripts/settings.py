### MARGARET GUO
### CREATED ON: 11/02/17 
### IMPORT ALL MODULES AND FILES NEEDED FOR the NCATS module

import sys, os
import pickle
# import csv
# import pandas as pd
import numpy as np
import time
# from scipy.stats import hypergeom
# import networkx as nx
# import numpy as np
# from collections import defaultdict
# import matplotlib
# import matplotlib.pyplot as plt
# import pandas as pd
# from textwrap import wrap
# from collections import defaultdict
# import itertools, time

# to import the rest of the modules more easily
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) # points to scripts DIR

from collections import defaultdict


# IMPORT USER_FUNCTIONS
from get_results import get_results
from get_output_data import get_output_data
import find_neighborhood_beta
from find_neighborhood_beta import find_neighborhood as fgn
import run_analysis as ann
from get_associations_deprecated import get_associations
import pharos_api
# import association file data
# Define variables with paths
global PARENT_DIR
global RSCS_DIR
global RESULTS_DIR
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # should yield the same as .. (with subfolders rcsc, scripts, and results)
RSCS_DIR = os.path.join(PARENT_DIR, 'rscs')
RESULTS_DIR = os.path.join(PARENT_DIR, 'results')


# loading data: interactome
NETWORKF = os.path.join(RSCS_DIR,'merged_interact_netx.pkl')#non-specific interaction network 

print('loading interactome data from', NETWORKF)
find_neighborhood_beta.GENE_GRAPH = pickle.load(open(NETWORKF,'rb'))
# print("(((((((((()))))))))")
# print(GENE_GRAPH.nodes())

# loading random data: 
ann.RAND_DIR = os.path.join(RESULTS_DIR, 'rand_iRefplus_intome/summary/')

# load drug targets, then format
DTF = os.path.join(RSCS_DIR, 'drug_intome_targets.pkl')
print('loading drug targets', DTF)
DTD = pickle.load(open(DTF,'rb'))


# load phenotype mapping data (map disease list of NCATs to possible phenotypes)
MAPPING_DICT = pickle.load(open(os.path.join(RSCS_DIR, 'conds_phen_matches_word_overlap.pkl') , 'rb'))
POSSIBLE_HEADER_RESULTS_FILE = ['phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']
# print('************')
# key = [(x,MAPPING_DICT[x]) for x in MAPPING_DICT.keys() if 'acne' in x]
# print(key,)
# print(len(key))
# Analysis parameters

SCORE_THRESHOLD_START = 0.8
SCORE_THRESHOLD_MIN = 0.5
SCORE_DELTA = 0.05