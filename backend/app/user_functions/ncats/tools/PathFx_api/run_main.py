'''
run_main.py

writen 7-13-17 JLW
Modified margaret, 10/26/17

Contains the relevant files necessary for an addendum function for Q2 --> having the ability to find drug indications (i.e. for drug repuposing),

'''
import sys, os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from settings import *

############################################ BASIC/SINGLE QUERIES #############################################

def preprocess_names(drug, disease=None):
    """
    does some form of preprocessing on drug and disease names so that they match the database

    inputs:
        <str> drug
        <str> disease (optional)
    outputs: tuple of
        <str> modified_drug_name
        <str> modified_disease_name (only if disease is defined)
    """
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
    DrugDisease.py calls this function:
    Given a drug, find all possible phenotypes it affects

        similar to run_drug_single except finds the significaent diseases of drug without referencing the given disease

    Input:
        <str> drug
        <str> storage_dir (optional)- where the resulting and intermediary analysis files are being stored
    Output
        <tuple> answer containing two strings
            ph_genes_str - a tab-delimited string with probability, Benjamin Hochberg significance cut off, the phenotype, and significant genes in subsequent separate columns
            drug - string of drug (return for debugging)

    '''
    #create storage directory
    results_storage_dir = os.path.join(RESULTS_DIR, storage_dir)
    if not os.path.exists(results_storage_dir):
        os.mkdir(results_storage_dir)

    drug = preprocess_names(drug)

    answer='None'
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
    outputs whether or not the disease was found in the drug network search 
    if yes, output signficant genes
    If no, iterate by lowering the threshold (default 0.8) until low point or found, 
        if still not found, output significant phenotypes found in the network and the genes associated with them

    Input
        <str> drug
        <str> disease
        <str> storage_dir
        <bool> iterate_until_found

    Output - a tuple of size 4 with 
        <bool> if the drug-disease relationship was found
        <str> explanation of found or not found, if not found a more categorical error msg
        <str> answer (see above)
        <str> drug (debug)
    '''
    results_storage_dir = os.path.join(RESULTS_DIR, storage_dir)
    if not os.path.exists(results_storage_dir):
        os.mkdir(results_storage_dir)


    drug, disease = preprocess_names(drug, disease)


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

            # POST PROCESSING 
            answer_found, answer = interpret_results(sum_asscs, disease)
            if  answer_found == 'found':
                print(answer, 'thres', threshold)

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


    else:
        print('Drug', drug, 'not in database')
        return False, ' '.join(['Drug', drug, 'not in database']), 'None', drug


################################### BATCH FUNCTIONS #################################


def run_drug_indications_multi(input_file, storage_dir = 'results_drugs_multi', save_file = 'results_drug_disease_matching.txt'):
    '''
    Run multi batches of find_drug_indications

    inputs:
        <str> input_file - path with drug (and disease) list
        <str>  storage_dir - directory where results are stored to
        <str>  save_file - the summary file name to be saved

    outputs:
        None - writes to files
    '''

    print('running drug_indications_multi from', input_file )
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
    answer_arr = []
    duration_arr = []
    for drug, disease in zip(drug_arr, disease_arr):
        start = time.time()
        answer, drug = find_drug_indications(drug, storage_dir)
        stop = time.time()
        duration = '%.2f' % (stop - start)
        answer_arr.append(answer)
        duration_arr.append(duration)
        print('>', drug, duration, answer)

    #get output as txt
    output_filename = os.path.join(RESULTS_DIR, storage_dir, save_file)
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'time_of_run(s)', 'answer\n')))

        for i in range(len(drug_arr)):
            line = '\t'.join((drug_arr[i], duration_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    print('written output to', output_filename)



def run_drug_multi(input_file, storage_dir = 'results_drugs_multi', save_file = 'results_drug_disease_matching.txt', iterate_until_found=False):
    '''
    Run multi batches of run_drug_single

    inputs:
        <str> input_file - path with drug (and disease) list
        <str>  storage_dir - directory where results are stored to
        <str>  save_file - the summary file name to be saved
        <bool> iterate_until_found - threshold lowering boolean

    outputs:
        None - writes to files
    '''

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
        duration = '%.2f' % (stop - start)
        answer_found_bool_arr.append(str(answer_found_bool))
        answer_found_arr.append(answer_found)
        answer_arr.append(answer)
        duration_arr.append(duration)
        print('>', drug, disease, duration, answer_found_bool, answer_found, answer)

    #get output as txt
    output_filename = os.path.join(RESULTS_DIR, storage_dir, save_file)
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'time_of_run(s)', 'answer_found_bool', 'answer_found', 'answer\n')))

        for i in range(len(drug_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], duration_arr[i], answer_found_bool_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    print('written output to', output_filename)


def interpret_results(sum_asscs, disease=None, mapping_dict = None):
    '''
    Takes in the information from the network analysis, formats the information in a form useful for output
    and if disease query is inputted, tries to find the phenotype that matches the query disease
    Input:
        <list> sum_accs: list of [ph,rank,BH,prb, asii,asin,pern,sig_genes] for each ph that is found by the algorithm to be significant plus the disease of interest
        <str> disease: being queried (optional), if None then just return all significant phenotypes 
        <dictionary> mapping_dict: dictionary of disease (keys) to list of possible phenotypes that the network produces (values) otherwise uses default dictionary in PathFx_DB
    Output:
        <str> interpretation of results 
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
            all_possible_diseases.append('\t'.join([ph, prb, BH ,', '.join(sig_genes)]))

    all_possible_diseases_string = '\t'.join(all_possible_diseases)
    if disease is None:
        return all_possible_diseases_string

    # check if there is a one-to-one match, returns three columns disease is the match we are looking for possible ph is found from network, the sig genes is list of genes
    # substrings are okay
    possible_diseases = [dis_key for dis_key in mapping_pkl.keys() if disease.lower() in dis_key]
    possible_phenotypes_matched = []
    for poss_dis in possible_diseases:
        possible_phenotypes_matched = possible_phenotypes_matched + mapping_pkl[poss_dis]

    if len(possible_phenotypes_matched ) > 0: # the disease is in the mapped list 
        # iterate through all posttibel of
        for possible_ph in possible_phenotypes_matched:
            # there was a phenotype match
            if possible_ph.lower().strip() in ph_to_sig_genes_dict.keys():
                sig_genes  = ph_to_sig_genes_dict[possible_ph.lower().strip()]
                print('\t'.join([prb, BH, disease, sig_genes]))
                return 'found', '\t'.join([prb, BH, disease, sig_genes])

        # if no match just return all possible phenotypes that it treats, where the first column is the probability 
        return ' '.join(['could not find disease - phenotype match in results']), all_possible_diseases_string
    else:
        return ' '.join(['could not find disease match in database']), all_possible_diseases_string




####################################### RETRIEVAL FUNCTIONS ##################################################

def retrieve_summary_per_drug_txt(drug, storage_dir):
    '''
    Finds information from a previous run in a storage directory

    Inputs:
        <str> drug
        <str> storage directory - hwere the relevant drug network query results are located
    Outputs:
        <sum_assc> - summary of associations to be fed into intepret_results

    '''

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

        return sum_asscs
    return None


def retrieve_results(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/ncats_test_intome_5_rerun_all_with_rand', old_format = False):
    """
    Reads in information from a run_drug_multi previous run to be interpretted by calling retrieve_summary_per_drug_txt

    Inputs:
        <str> input_file - query file with tab-delimited list of drug-disease pairs
        <str> storage_dir - where previous run to search through is 
        <bool> old_format - based on a version change in the original code

    """
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
    
    # iterate through given mapping tdict
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
    with open(output_filename, 'w')  as f:
        f.write('\t'.join(('drug', 'condition', 'answer_found', 'answer\n')))

        for i in range(len(answer_found_arr)):
            line = '\t'.join((drug_arr[i], disease_arr[i], answer_found_arr[i], answer_arr[i]))
            f.write(''.join((line, '\n')))

    
    print('written output to', output_filename)


# unit testing:
if __name__ == '__main__':

    retrieve_results(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/ncats_test_intome_5_rerun_all_with_rand', old_format = True)
    run_drug_single('Latamoxef', 'Bacterial Infections', storage_dir = "../results/results_test")
    run_drug_multi(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/results_drugs_multi_time2', save_file = 'results_drug_disease_matching.txt', iterate_until_found=False)
    run_drug_indications_multi(input_file = '../rscs/q2-drugandcondition-list.txt', storage_dir = '../results/results_drugs_indications_multi_time', save_file = 'results_drug_disease_matching.txt')
    analysis_name = 'ncats_test_intome_6_remove24_min_networkassoc' 
    rdir = '../results/'+analysis_name+'/'
    if not os.path.exists(rdir):
        os.mkdir(rdir)

