### MARGARET GUO
### CREATED ON: 11/02/17 
### IMPORT ALL MODULES AND FILES NEEDED FOR the run_main.py

import sys, os
import pickle
import numpy as np
import time


# to import the rest of the modules more easily
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) # points to scripts DIR

from collections import defaultdict


# IMPORT USER_FUNCTIONS
from app.user_functions.ncats.tools.Pharos_api import pharos_api

from get_results import get_results
from get_output_data import get_output_data
import find_neighborhood_beta
from find_neighborhood_beta import find_neighborhood as fgn
import run_analysis as ann
from get_associations import get_associations


# import association file data
# Define variables with paths
global PARENT_DIR
global RSCS_DIR
global RESULTS_DIR
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # should get to the ncats directory
RSCS_DIR = os.path.join(PARENT_DIR, 'DB_data/PathFx_DB')
RESULTS_DIR = os.path.join(PARENT_DIR, 'results')


# loading data: interactome
NETWORKF = os.path.join(RSCS_DIR,'merged_interact_netx.pkl')#non-specific interaction network 
print('loading interactome data from', NETWORKF)
find_neighborhood_beta.GENE_GRAPH = pickle.load(open(NETWORKF,'rb'))


# loading random data: 
ann.RAND_DIR = os.path.join(RSCS_DIR, 'rand_iRefplus_intome/summary/')

# load drug targets, then format
DTF = os.path.join(RSCS_DIR, 'drug_intome_targets.pkl')
print('loading drug targets', DTF)
DTD = pickle.load(open(DTF,'rb'))


# load phenotype mapping data (map disease list of NCATs to possible phenotypes)
MAPPING_DICT = pickle.load(open(os.path.join(RSCS_DIR, 'conds_phen_matches_word_overlap.pkl') , 'rb'))
POSSIBLE_HEADER_RESULTS_FILE = ['phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']


SCORE_THRESHOLD_START = 0.8
SCORE_THRESHOLD_MIN = 0.5
SCORE_DELTA = 0.05