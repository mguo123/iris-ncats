# written to create a better summary table
# from phenotypic enrichment, written
# 8-31-17 JLW
# Modified margaret, 10/26/17


import csv, os, pickle
import numpy as np
from collections import defaultdict

# Setting paths
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # should get to the ncats directory
RSCS_DIR = os.path.join(PARENT_DIR, 'DB_data/PathFx_DB')
MAPPING_DICT = pickle.load(open(os.path.join(RSCS_DIR, 'conds_phen_matches_word_overlap.pkl') , 'rb'))
POSSIBLE_HEADER_RESULTS_FILE = ['phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']


def get_sum(mf,gtphen, disease=None): 
    """
    Pass a merged phenotype result file and get summary of results 

    Inputs:
        <str> mf: merged phenotyp result file path
        <gtphen> path to pkl of phenotype to gene information
        <str> disease: optional
    
    Outputs:
        <list> sig_phens_data - list of list of strings into information to be printed to file in get_results
    """

    # pass a merged phenotype result file
    sig_phens_data = []
    gtpd = pickle.load(open(gtphen,'rb'))
    dR = csv.DictReader(open(mf,'rU', encoding='utf-8'),delimiter='\t', quotechar='"')
    print(dR)


    possible_phenotypes_matched = []
    if disease is not None:
        possible_diseases = [dis_key for dis_key in MAPPING_DICT.keys() if disease.lower() in dis_key]
        for poss_dis in possible_diseases:
            possible_phenotypes_matched = possible_phenotypes_matched + MAPPING_DICT[poss_dis]

    for l in dR:
        # gather relevant data
        [ph,rank,BH,prb] = [l['phenotype'],l['rank'],l['Benjamini-Hochberg'],l['probability']]
        [asii,asin] = [l['assoc in intom'],l['assoc in neigh']]    
        pern = float(asin)/float(asii)

        # determine significance
        if float(prb) < float(BH) and float(asii)>10:
            sig_genes = gtpd[ph]
            sig_phens_data.append((ph,rank,BH,prb, asii,asin,str(pern),sig_genes))
        elif len(possible_phenotypes_matched) > 0: # the disease is in the mapped list 
            if ph.lower().strip() in possible_phenotypes_matched: # if the phenotype corresponds to the disease in question
                sig_genes = gtpd[ph]
                sig_phens_data.append((ph,rank,BH,prb, asii,asin,str(pern),sig_genes))
                break # exit loop as soon as find the phenotype


    return sig_phens_data


def get_results(storage_dir, drug, disease=None, old_format = False, save_file = True):
    """
    Get phenotype information from the merged neighborhood associations and writes to a tab-delimited file containing:
    phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen

    Inputs:
        <str> storage_dir
        <str> drug
        <str> disease (optional)
        <bool> old_format - just use defaults, for debugging
        <bool> save_file - for debugging, default True

    Outputs:
        <list> sum_asscs - significant phenotype information stored as a list of of lists containing
            (all <str>)
            ph: phenotype
            rank: rank
            BH: BH threshold value
            prb: p-value, 
            acii: association in intom
            acin: association in neighborhood
            pern: percent overlap of neighboorhoods between drug and phenotype
            sig_genes: number of neighboring genes in phen
    """

    # getting analysis info from file
    drug_dir = os.path.join(storage_dir, drug + "_networks")
    if old_format:
        drug_file = os.path.join(drug_dir, drug + '_merged_neighborhood_assoc_prb.txt')

    else:
        drug_file = os.path.join(drug_dir, drug + '_merged_neighborhood_assocations.txt')
    print(drug)

    # if you can find all the relevant inforamtion 
    if os.path.isfile(drug_file): 

        # get phenotype signficance information
        mapf = drug_file
        gtphen = os.path.join(drug_dir,drug+'_phens_to_genes.pkl')
        sum_asscs = get_sum(mapf,gtphen, disease)


        # write to file
        outfname = os.path.join(drug_dir,drug+'_merged_assc_full_summary.txt')
        
        if save_file:
            outf = open(outfname,'w')

            print('writing output to', outfname)
            outf.write('\t'.join(['phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']))

            for line in sum_asscs:
                if line is not None:
                    [ph,rank,BH,prb, asii,asin,pern,sig_genes] = line
                    sig_gn_str = ','.join(sig_genes)
                    outf.write('\t'.join([ph,rank,BH,prb, asii,asin,pern,sig_gn_str,'\n']))        
            outf.close()

        return sum_asscs

    else:
        return []


# if __name__ == "__main__":
#         main()
