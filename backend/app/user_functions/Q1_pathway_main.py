

import subprocess
import os

import pandas as pd

# sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) # points to E DIR

overall_path = os.path.abspath(os.path.dirname(__file__))
results_dir = os.path.join(overall_path, "ncats/results/Q2/")

if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# from app.user_functions.ncats.EBC_api import EBC_api
# from app.user_functions.ncats.Pharos_api import pharos_api
# from app.user_functions.ncats.GO_api import go_api

from ncats.EBC_api import EBC_api
from ncats.Pharos_api import pharos_api
from ncats.GO_api import go_api

"""
Q2_query

Runs the Q2 pipeline with the given drug and disease

Input:
QDrug - String of the drug common name
QDisease - String of the disease common name
gen_image - Whether you would like the script to render a network graphic

Output:
results - Pandas dataframe of GO Enrichment for the drug-disease combo with the following columns:
'M'
N' 
'n' 
'name' 
'namespace'
'p'
'q'
'rejected'
'term'
'x'

"""

def Q1_pathway_query(QDisease1, QDisease2, gen_tissues_image=False, gen_interaction_image=False, output_full=False):

    # Pre-process query
    disease1 = QDisease1.strip().lower()

    disease2 = QDisease2.strip().lower()

    # Generate prefix for output file, is drug_disease
    out_name = disease1.replace(" ", "-").lower() + "_" + disease2.replace(" ", "-").lower()

    # Get list of genes causally annotated to a disease
    dis_gene_list1 = EBC_api.get_disease_gene_list(disease1, freq_correct=True)
    dis_gene_list2 = EBC_api.get_disease_gene_list(disease2, freq_correct=True)
    # print(dis_gene_list1, dis_gene_list2)
    # # If first query did not work, try getting nearest matches for the query and try again
    if dis_gene_list1 is None:

        # Get matches
        d_temp, match_type = EBC_api.query_term_for_matches(disease1)

        # If not fuzzy matching, but simple order matching, use the new disease query and proceed
        if match_type == 0:
            dis_gene_list1 = EBC_api.get_disease_gene_list(d_temp, freq_correct=True)

    # If first query did not work, try getting nearest matches for the query and try again
    if dis_gene_list2 is None:

        # Get matches
        d_temp, match_type = EBC_api.query_term_for_matches(disease2)

        # If not fuzzy matching, but simple order matching, use the new disease query and proceed
        if match_type == 0:
            dis_gene_list2 = EBC_api.get_disease_gene_list(d_temp, freq_correct=True)

    # If either disease or drug list comes up empty
    # the query has failed and we return a statement to that effect

    if dis_gene_list1 is None and dis_gene_list2 is None:
        return "Diseases not recognized, better luck next time..."

    elif dis_gene_list1 is None:
        return "Disease 1 not recognized"

    elif dis_gene_list2 is None:
        return "Disease 2 not recognized"

    # If we have targets
    else:

        #############################
        # PMIDs = EBC_api.query_chemical_disease(drug, disease, get_PMIDs=True)
        # print(PMIDs)
        #############################


        # Select the top 25 genes from the disease gene list for GO enrichment
        dis_genes1 = [[EBC_api.resolve_EntrezGeneID_to_NCBIGeneName(x),x] for x in dis_gene_list1]

        dis_genes1 = pd.DataFrame(dis_genes1, columns=["Gene", "Entrez ID"])
        dis_gene_list1 = list(map(int, dis_gene_list1))

        # print(dis_genes1)

        dis_genes2 = [[EBC_api.resolve_EntrezGeneID_to_NCBIGeneName(x),x] for x in dis_gene_list2]

        dis_genes2 = pd.DataFrame(dis_genes2, columns=["Gene", "Entrez ID"])
        dis_gene_list2 = list(map(int, dis_gene_list2))

        # print(dis_genes2)

        # Get GO Enrichment statistics
        result1 = GO_API.calc_GO_enrichment(dis_gene_list1, os.path.join(results_dir, out_name))
        result2 = GO_API.calc_GO_enrichment(dis_gene_list2, os.path.join(results_dir, out_name))


        if output_full is False:
            result1 = result1.loc[result1['rejected'] == 1.0, ['namespace', 'name', 'p', 'q', 'gene_target']]
            result1 = result1.sort_values(by=['q'])
            result2 = result2.loc[result2['rejected'] == 1.0, ['namespace', 'name', 'p', 'q', 'gene_target']]
            result2 = result2.sort_values(by=['q'])

        A = set(list(result1['name']))
        # print(result1.shape)

        B = set(list(result2['name']))
        # print(result2.shape)
        print(B)
        print(len(A), len(B))
        print(len(A.intersection(B)))
        print(float(len(A.intersection(B))) / float(len(A) + len(B)))

            #
        if gen_interaction_image:
            subprocess.check_call(['dot', '-Tpng', os.path.join(results_dir, out_name)+ '.dot', '-o', os.path.join(results_dir, out_name) + '.png'])
        if gen_tissues_image:
            pass


        # result = {"GOENRICH":result, "drug_genes":drug_genes, "disease_genes":dis_genes}

    # return(drug_tissues)
    # return(result)

# unit testing
if __name__ == "__main__":
    disease1 = "malaria"

    disease2 = "sickle cell anemia"

    GO_API = go_api.GO_api("./ncats/GO_api/GO_DB")

    print(disease1,disease2)
    # resultlt = \
    Q1_pathway_query(disease1, disease2)


    disease1 = "malaria"

    disease2 = "breast cancer"

    print(disease1,disease2)
    # resultlt = \
    Q1_pathway_query(disease1, disease2)

    # if type(result) is not str:
    #     print(result)
    # else:
    #     print(result)

