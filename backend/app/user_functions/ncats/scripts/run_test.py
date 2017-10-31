    # written as a template to run a class of drugs and create a heatmap of results
# writen 7-13-17 JLW
# Modified margaret, 10/26/17


# Analysis parameters
#### FIX ####
import csv
import pandas as pd
import numpy as np

SCORE_THRESHOLD_START = 0.8
SCORE_THRESHOLD_MIN = 0.5
SCORE_DELTA = 0.05



def run_drug_single(drug, disease, storage_dir = 'results_drugs', iterate_until_found=False):
    return 'False', 'None'


def run_drug_multi(input_file, storage_dir = 'results_drugs_multi', iterate_until_found=False):
    print('running drug_multi from', input_file )
    # get in input_file of two columns of drugs and diseases, that shouldbe paired together
    if input_file.endswith('.txt'):
        with open(input_file, 'r') as f:
            x= f.readlines()
            x = np.array([line.strip().split('\t') for line in x][1:])
            drug_arr = x[:,0]
            disease_arr = x[:,1]
    elif input_file.endswith('.csv'):
        data = pd.read_csv(input_file, sep=',', quotechar='\"')
        drug_arr = np.array(data[list(data)[0]])
        disease_arr = np.array(data[list(data)[1]])
    print(drug_arr)
    print(disease_arr)

run_drug_multi('/Users/margaret/Documents/iris-agent/node_modules/ncats/rscs/q2-drugandcondition-list.csv', storage_dir = 'abc', iterate_until_found=False)