    # written as a template to run a class of drugs and create a heatmap of results
# writen 7-13-17 JLW
# Modified margaret, 10/26/17


import csv, pickle, os, sys, find_neighborhood_beta
import networkx as nx
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import pandas as pd
from textwrap import wrap

# from optparse import OptionParser
from collections import defaultdict

from get_results import get_results
from get_output_data import get_output_data
from find_neighborhood_beta import find_neighborhood as fgn
import run_analysis as ann

# Analysis parameters
#### FIX ####

SCORE_THRESHOLD_START = 0.8
SCORE_THRESHOLD_MIN = 0.5
SCORE_DELTA = 0.05

def preprocess_names(drug, disease):

    drug = drug.lower()
    first_letter = drug[0].upper()
    modified_drug_name = first_letter + drug[1:]
    modified_disease = disease.lower().strip('\"')
    return modified_drug_name, modified_disease

def run_drug_single(drug, disease, storage_dir = 'results_drugs', iterate_until_found=False):
    if not os.path.exists(storage_dir):
        os.mkdir(storage_dir)


    drug, disease = preprocess_names(drug, disease)

    #### TODO: make global list of paths
    # loading data: interactome
    networkf = '../rscs/merged_interact_netx.pkl' #non-specific interaction network #### CHANGE LOCATION
    print('loading interactome data', networkf)
    GENE_GRAPH = pickle.load(open(networkf,'rb'))

    # loading random data: 
    ann.rand_dir = '../results/rand_iRefplus_intome/summary/'

    # load drug targets, then format
    dtf = '../rscs/drug_intome_targets.pkl'
    print('loading drug targets', dtf)
    dtd = pickle.load(open(dtf,'rb'))
    print('running analysis of', drug, 'to treat', disease)
    if drug in dtd:
        dts = [(drug.replace(' ',''), dtd[drug])] # remove white spaces for later analysis

        # create networks, merge, do phenotype enrichment
        threshold = SCORE_THRESHOLD_START # STARTING THRESHOLD
        answer_found = False
        while not answer_found and threshold >= SCORE_THRESHOLD_MIN :
            all_merge_files = ann.run_all_drugs(dts,storage_dir,GENE_GRAPH,threshold)


            # get/store output data
            assoc_to_genes = get_output_data(storage_dir, drug)
            sum_asscs = get_results(storage_dir, drug)
            # print('here!!!', sum_asscs)

    #             #### POST PROCESSING ####
            answer_found, answer = interpret_results(sum_asscs, disease)
            if answer_found:
                print(answer, 'thres', threshold)
                ' '.join((str(threshold), answer))
                #'phenotype','rank','BHcorrPval','assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen', 'threshold'
                return str(answer_found), answer 
            else:
                if iterate_until_found:
                    threshold -= SCORE_DELTA
                    print('new threshold', threshold)
                else:
                    threshold = -1

        print( disease, 'not found to be treated by', drug, 'but found possible matches', answer)
        return str(answer_found), answer
#     z = pickle.load(open('conds_phen_matches_word_overlap.pkl', 'rb'))
# for key, value in z.items():
#     print(key, value)

    else:
        print('Drug', drug, 'not in database')
        return 'False', 'None'


def run_drug_multi(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = 'results_drugs_multi', iterate_until_found=False):
    print('running drug_multi from', input_file )
    # get in input_file of two columns of drugs and diseases, that shouldbe paired together
    if input_file.endswith('.txt'):
        with open(input_file, 'r') as f:
            x= f.readlines()
            x = np.array([line.strip().split('\t') for line in x][1:])
            drug_arr = x[:,0]
            disease_arr = x[:,1]
    elif input_file.endswith('.csv'):
        data = pd.read_csv(input_file, sep=',')
        drug_arr = np.array(data[list(data)[0]])
        disease_arr = np.array(data[list(data)[1]])

    #iterate through list
    answer_found_arr = []
    answer_arr = []
    for drug, disease in zip(drug_arr, disease_arr):
        answer_found, answer = run_drug_single(drug, disease, storage_dir, iterate_until_found)
        answer_found_arr.append(answer_found)
        answer_arr.append(answer)


    #get output as txt
    output_filename = os.path.join(storage_dir, 'results_drug_disease_matching.txt')
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'answer_found', 'answer\n')))

        for i in range(len(drug_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    print('written output to', output_filename)


def interpret_results(sum_asscs, disease, mapping_dict = None):
    if mapping_dict is None:
        mapping_pkl = '../rscs/conds_phen_matches_word_overlap.pkl'
        mapping_dict = pickle.load(open(mapping_pkl, 'rb'))

    answer_found = False


    ph_to_sig_genes_dict = defaultdict(list)
    all_possible_diseases = []
    for line in sum_asscs:
        if len(line) == 7:
            [ph,rank,BH,asii,asin,pern,sig_genes] = line
            ph_to_sig_genes_dict[ph.lower().strip()] = ''.join(sig_genes)
            all_possible_diseases.append('\t'.join([BH, ph,sig_genes]))
 
    # check if there is a one-to-one match, returns three columns disease is the match we are looking for possible ph is found from network, the sig genes is list of genes
    possible_phenotypes_matched = mapping_dict[disease.lower()]
    # print('possible_phenotypes_matched  ', possible_phenotypes_matched  )
    for possible_ph in possible_phenotypes_matched:
        if possible_ph.lower().strip() in ph_to_sig_genes_dict.keys():
            answer_found = True
            sig_genes  = ph_to_sig_genes_dict[possible_ph.lower().strip()]
            print('\t'.join([disease, possible_ph, sig_genes]))
            return answer_found, '\t'.join([disease, possible_ph, sig_genes])
    
    # if no match just return all possible phenotypes that it treats, where the first column is the probability 
    all_possible_diseases_string = '\t'.join(all_possible_diseases)
    return answer_found, all_possible_diseases_string


# def reverse_mapping_dict(mapping_dict):
#     # turn disease to possible phenotypes dictionary to a phenotype to possible disease
#     ph_to_disease_dict = defaultdict(list)
#     for disease, ph_list in mapping_dict.items():
#         for ph in ph_list:
#             ph_to_disease_dict[ph].append(disease)

## for getting from an actual folder
def retrieve_summary_per_drug_txt(drug, storage_dir):

    outfname = os.path.join(storage_dir, drug + "_networks",drug+'_merged_assc_full_summary.txt')
    # print(outfname)
    if os.path.isfile(outfname):
        with open(outfname,'r') as f:
            sum_asscs = f.readlines()
        # sum_asscs = []
        # print(lines)
        for i, line in enumerate(sum_asscs):
            sum_asscs[i] = line.strip().split('\t')

        if len(sum_asscs ) > 0:
            sum_asscs.pop(0)
            # print(line)
            # sum_asscs.append(line)

        # outf.red('\t'.join(['phenotype','rank','BHcorrPval','assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']))
        # for [ph,rank,BH,asii,asin,pern,sig_genes] in sum_asscs:
        #     sig_gn_str = ','.join(sig_genes)
        #     outf.write('\t'.join([ph,rank,BH,asii,asin,pern,sig_gn_str,'\n']))        
        # outf.close()
        return sum_asscs
    return None

def retrieve_results(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/ncats_test_intome_5_rerun_all_with_rand', old_format = False):
    print('retrieving results from', storage_dir)
    # get in input_file of two columns of drugs and diseases, that shouldbe paired together
    with open(input_file, 'r') as f:
        x= f.readlines()
    drug_arr = []
    disease_arr = []
    for i, line in enumerate(x[1:]): # skip header
        drug, disease = line.strip().split('\t')
        drug, disease = preprocess_names(drug, disease)
        drug_arr.append(drug)
        disease_arr.append(disease)
    
    # Get mapping tdict
    mapping_pkl = '../rscs/conds_phen_matches_word_overlap.pkl'
    mapping_dict = pickle.load(open(mapping_pkl, 'rb'))
    #iterate through list
    answer_found_arr = []
    answer_arr = []
    for i, (drug, disease) in enumerate(zip(drug_arr, disease_arr)):
        # if i == 65:
        if i < len(drug_arr):
            print('getting results', i, drug, disease)
            sum_asscs = retrieve_summary_per_drug_txt(drug, storage_dir)
            # sum_asscs = get_results(storage_dir, drug, old_format = old_format,  save_file = False) #True)   
            if sum_asscs is not None: 
                answer_found, answer = interpret_results(sum_asscs, disease,mapping_dict)
            else:

                answer_found, answer = 'False', 'No diseases detected in network'

            answer_found_arr.append(str(answer_found))
            answer_arr.append(answer)
        
    #get output as txt
    output_filename = os.path.join(storage_dir, 'results_drug_disease_matching_test.txt')
    # output_filename = os.path.join(storage_dir, 'results_drug_disease_matching.txt')
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'answer_found', 'answer\n')))

        for i in range(len(answer_found_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    
    print('written output to', output_filename)

##### STUFF I'm RUNNING
# retrieve_results(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/ncats_test_intome_5_rerun_all_with_rand', old_format = True)
# run_drug_multi(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = 'results_drugs_multi2', iterate_until_found=False)

# analysis_name = 'ncats_test_intome_6_remove24_min_networkassoc' 
# rdir = '../results/'+analysis_name+'/'
# if not os.path.exists(rdir):
#     os.mkdir(rdir)


#### SAVE FOR MULTI_DRUG and disease inputs, right now for interactive purposes should do one by one.
# ####add or load drug list from NCATs #### 
# with open('../rscs/q2-drugandcondition-list.txt', 'r') as f:
#     x= f.readlines()
# drugs = []
# diseases = []
# for i, line in enumerate(x[1:]):
#     line = line.strip()
#     drug, disease = line.split('\t')
#     drug = drug.lower()
#     first_letter = drug[0].upper()
#     modified_drug_name = first_letter + drug[1:]
#     drugs.append(modified_drug_name)
#     diseases.append(disease)


###### ADD in drug names that Jenn hasn't modelled
#with open('../rscs/drug_list_not_modelled.txt', 'r') as f:
#    x= f.readlines()
#drugs = []
#for i, line in enumerate(x[1:]):
#    drug = line.strip()
#    drug = drug.lower()
#    first_letter = drug[0].upper()
#    modified_drug_name = first_letter + drug[1:]
#    drugs.append(modified_drug_name)
#
# print(drugs)
# dtf = '../rscs/drug_intome_targets.pkl'
# dtd = pickle.load(open(dtf,'rb'))
# print(dtd)
# dts = [(d,dtd[d]) for d in drugs if d in dtd]
# print(dts)
# dts = [(d.replace(' ',''),tlist) for (d,tlist) in dts] # remove white spaces for later analysis


### IGNORE ME #### check to see if all the drugs have targets in the network
# print(len(dts), 'len dts')
# ann.check_if_drug_in_network(dts, GENE_GRAPH)
# raise ValueError("TESTING")
# create networks, merge, do phenotype enrichment
# all_merge_files = ann.run_all_drugs(dts,rdir,GENE_GRAPH,SCORE_THRESHOLD)



### IGNORED
# call method to plot a heatmap
# per_cutoff = 0.05 #can change this to include more drugs in the heatmap 
# ann.make_a_heatmap(drugs,all_merge_files,rdir,analysis_name,per_cutoff)

        
# #### HACK FOR NCATS
# os.system(' '.join(('python get_genes_to_phen_dics.py -d', rdir)))
# os.system(' '.join(('python summarize_network_phenotypes.py -d', rdir)))

