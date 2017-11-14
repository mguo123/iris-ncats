    # written as a template to run a class of drugs and create a heatmap of results
# writen 7-13-17 JLW
# Modified margaret, 10/26/17
import sys, os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from settings import *



# import csv, pickle, os, sys
# sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
# print(sys.path)
# import networkx as nx
# import numpy as np
# from collections import defaultdict
# import matplotlib
# # matplotlib.use("AGG")
# import matplotlib.pyplot as plt
# import pandas as pd
# from textwrap import wrap

# # from optparse import OptionParser
# from collections import defaultdict

# from get_results import get_results
# from get_output_data import get_output_data
# import find_neighborhood_beta
# from find_neighborhood_beta import find_neighborhood as fgn
# import run_analysis as ann

# # Analysis parameters
# #### FIX ####

# SCORE_THRESHOLD_START = 0.8
# SCORE_THRESHOLD_MIN = 0.5
# SCORE_DELTA = 0.05

def preprocess_names(drug, disease=None):

    drug = drug.lower()
    first_letter = drug[0].upper()
    modified_drug_name = first_letter + drug[1:]
    if disease is not None:
        modified_disease = disease.lower().strip('\"')
        return modified_drug_name, modified_disease
    else:
        return modified_drug_name

def find_drug_indications(drug, storage_dir = 'Q2'):
    '''
    similar to run_drug single except finds the significaent diseases of drug
    '''
    #create storage directory
    results_storage_dir = os.path.join(RESULTS_DIR, storage_dir)
    if not os.path.exists(results_storage_dir):
        os.mkdir(results_storage_dir)

    drug = preprocess_names(drug)

    if drug in DTD:
        dts = [(drug.replace(' ',''), DTD[drug])] # remove white spaces for later analysis

        # create networks, merge, do phenotype enrichment
        threshold = SCORE_THRESHOLD_START # STARTING THRESHOLD

        # get the merged files
        all_merge_files = ann.run_all_drugs(dts,results_storage_dir,find_neighborhood_beta.GENE_GRAPH,threshold)


        # get/store output data
        assoc_to_genes = get_output_data(results_storage_dir, drug)
        sum_asscs = get_results(results_storage_dir, drug)
        # print('here!!!', sum_asscs)

#             #### POST PROCESSING ####
        answer = interpret_results(sum_asscs)
        
        return answer, drug

    else:
        return answer, drug


def run_drug_single(drug, disease, storage_dir = 'Q2', iterate_until_found=False):
    '''
    Takes in strings for drug and disease then 
    outputs the phenotypes found in the network and the genes associated with them
    '''
    results_storage_dir = os.path.join(RESULTS_DIR, storage_dir)
    if not os.path.exists(results_storage_dir):
        os.mkdir(results_storage_dir)


    drug, disease = preprocess_names(drug, disease)

    ##### ALL OF THESE ARE DEFINED IN settings.py as global vars
    # # loading data: interactome
    # NETWORKF = '../rscs/merged_interact_netx.pkl' #non-specific interaction network #### CHANGE LOCATION
    # print('loading interactome data', NETWORKF)
    # GENE_GRAPH = pickle.load(open(NETWORKF,'rb'))

    # # loading random data: 
    # ann.rand_dir = '../results/rand_iRefplus_intome/summary/'

    # # load drug targets, then format
    # DTF = '../rscs/drug_intome_targets.pkl'
    # print('loading drug targets', DTF)
    # DTD = pickle.load(open(DTF,'rb'))

    print('running analysis of', drug, 'to treat', disease)
    if drug in DTD:
        dts = [(drug.replace(' ',''), DTD[drug])] # remove white spaces for later analysis

        # create networks, merge, do phenotype enrichment
        threshold = SCORE_THRESHOLD_START # STARTING THRESHOLD
        answer_found = False
        while not answer_found and threshold >= SCORE_THRESHOLD_MIN :
            all_merge_files = ann.run_all_drugs(dts,results_storage_dir,find_neighborhood_beta.GENE_GRAPH,threshold)


            # get/store output data
            assoc_to_genes = get_output_data(results_storage_dir, drug)
            sum_asscs = get_results(results_storage_dir, drug, disease)
            # print('here!!!', sum_asscs)

    #             #### POST PROCESSING ####
            answer_found, answer = interpret_results(sum_asscs, disease)
            if  answer_found == 'found':
                print(answer, 'thres', threshold)
                ' '.join((str(threshold), answer))
                #'phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen', 'threshold'
                return True, answer_found, answer, drug
            else:
                if iterate_until_found:
                    threshold -= SCORE_DELTA
                    print('new threshold', threshold)
                else:
                    threshold = -1

        print( disease, 'not found to be treated by', drug, 'but found possible matches')
        return False, str(answer_found), answer, drug
#     z = pickle.load(open('conds_phen_matches_word_overlap.pkl', 'rb'))
# for key, value in z.items():
#     print(key, value)

    else:
        print('Drug', drug, 'not in database')
        return False, ' '.join(['Drug', drug, 'not in database']), 'None', drug


def run_drug_multi(input_file, storage_dir = 'results_drugs_multi', save_file = 'results_drug_disease_matching.txt', iterate_until_found=False):
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
    answer_found_bool_arr = []
    answer_found_arr = []
    answer_arr = []
    duration_arr = []
    for drug, disease in zip(drug_arr, disease_arr):
        start = time.time()
        answer_found_bool, answer_found, answer, drug = run_drug_single(drug, disease, storage_dir, iterate_until_found)
        stop = time.time()
        duration = '%.1f' % (stop - start)
        answer_found_bool_arr.append(answer_found_bool)
        answer_found_arr.append(answer_found)
        answer_arr.append(answer)
        duration_arr.append(duration)


    #get output as txt
    output_filename = os.path.join(RESULTS_DIR, storage_dir, save_file)
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'time_of_run(s)', 'answer_found_bool', 'answer_found', 'answer\n')))

        for i in range(len(drug_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], duration[i], answer_found_bool_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    print('written output to', output_filename)


def interpret_results(sum_asscs, disease=None, mapping_dict = None):
    '''
    sum_accs: list of [ph,rank,BH,prb, asii,asin,pern,sig_genes] for each ph that is found by the algorithm to be significant plus the disease of interest
    mapping_dict: dictionary of disease (keys) to list of possible phenotypes that the network produces (values)
    RETURNS!!!
    '''
    print('interpretting results')
    if mapping_dict is None:
        mapping_pkl = MAPPING_DICT


    ph_to_sig_genes_dict = defaultdict(list)
    all_possible_diseases = []
    for line in sum_asscs:
        if len(line) == len(POSSIBLE_HEADER_RESULTS_FILE):
            [ph,rank,BH,prb, asii,asin,pern,sig_genes] = line
            BH = '%.3g' % float(BH)
            prb = '%.3g' % float(prb)
            ph_to_sig_genes_dict[ph.lower().strip()] = ','.join(sig_genes)
            all_possible_diseases.append('\t'.join([prb, BH, ph,', '.join(sig_genes)]))

    all_possible_diseases_string = '\t'.join(all_possible_diseases)
    if disease is None:
        return all_possible_diseases_string

    # check if there is a one-to-one match, returns three columns disease is the match we are looking for possible ph is found from network, the sig genes is list of genes
    # substrings are okay
    possible_diseases = [dis_key for dis_key in MAPPING_DICT.keys() if disease.lower() in dis_key]
    possible_phenotypes_matched = []
    for poss_dis in possible_diseases:
        possible_phenotypes_matched = possible_phenotypes_matched + MAPPING_DICT[poss_dis]

    if len(possible_phenotypes_matched ) > 0: # the disease is in the mapped list 
        for possible_ph in possible_phenotypes_matched:
            if possible_ph.lower().strip() in ph_to_sig_genes_dict.keys():
                sig_genes  = ph_to_sig_genes_dict[possible_ph.lower().strip()]
                # print(possible_ph, sig_genes, 'FOUND IN possible_phenotypes_matched')
                print('\t'.join([prb, BH, disease, sig_genes]))
                return 'found', '\t'.join([prb, BH, disease, sig_genes])
            # print(possible_ph)

        # if no match just return all possible phenotypes that it treats, where the first column is the probability 
        
        return ' '.join(['could not find disease - phenotype match in results']), all_possible_diseases_string
    else:
        return ' '.join(['could not find disease match in database']), all_possible_diseases_string


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
    # mapping_pkl = '../rscs/conds_phen_matches_word_overlap.pkl'
    # mapping_dict = pickle.load(open(mapping_pkl, 'rb'))
    #iterate through list
    answer_found_arr = []
    answer_arr = []
    for i, (drug, disease) in enumerate(zip(drug_arr, disease_arr)):
        # if i == 47:  #### DEBUG
        if i < len(drug_arr):
            print('getting results', i, drug, disease)
            # sum_asscs = retrieve_summary_per_drug_txt(drug, storage_dir)
            sum_asscs = get_results(storage_dir, drug, old_format = old_format,  save_file = False) #True)   
            if sum_asscs is not None: 
                answer_found, answer = interpret_results(sum_asscs, disease,MAPPING_DICT)
            else:

                answer_found, answer = 'False', 'No diseases detected in network'

            answer_found_arr.append(str(answer_found))
            answer_arr.append(answer)
        
    #get output as txt
    output_filename = os.path.join(storage_dir, 'results_drug_disease_matching_test2.txt')
    # output_filename = os.path.join(storage_dir, 'results_drug_disease_matching.txt')
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'answer_found', 'answer\n')))

        for i in range(len(answer_found_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    
    print('written output to', output_filename)

##### STUFF I'm RUNNING
# retrieve_results(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/ncats_test_intome_5_rerun_all_with_rand', old_format = True)
# run_drug_multi(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/results_drugs_multi2', save_file = 'results_drug_disease_matching.txt', iterate_until_found=False)

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


