

import subprocess
import os

import pandas as pd

# sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) # points to E DIR

overall_path = os.path.abspath(os.path.dirname(__file__))
results_dir = os.path.join(overall_path, "ncats/results/Q2/")

if not os.path.exists(results_dir):
    os.makedirs(results_dir)

from app.user_functions.ncats.EBC_api import EBC_api
from app.user_functions.ncats.Pharos_api import pharos_api
from app.user_functions.ncats.GO_api import go_api

GO_API = go_api.GO_api(os.path.join(overall_path, "ncats/GO_api/GO_DB"))

# from ncats.EBC_api import EBC_api
# from ncats.Pharos_api import pharos_api
# from ncats.GO_api import go_api

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

def Q2_query(QDrug, QDisease, gen_tissues_image=False, gen_interaction_image=False, output_full=False):
    #
    # Pre-process query
    drug = QDrug.strip().lower()

    disease = QDisease.strip().lower()

    # Generate prefix for output file, is drug_disease
    out_name = drug + "_" + disease.replace(" ", "-").lower()

    # Get list of genes causally annotated to a disease
    dis_gene_list = EBC_api.get_disease_gene_list(disease, freq_correct=True)



    # If first query did not work, try getting nearest matches for the query and try again
    if dis_gene_list is None:

        # Get matches
        disease2, match_type = EBC_api.query_term_for_matches(disease)

        # If not fuzzy matching, but simple order matching, use the new disease query and proceed
        if match_type == 0:
            dis_gene_list = EBC_api.get_disease_gene_list(disease2, freq_correct=True)

    # Get list of drug targets from Pharos
    drug_genes = pharos_api.get_ligand_targets(drug)
    # # Get list of tissues from Pharos
    # drug_tissues = pharos_api.get_drug_tissues(drug)

    # If Pharos did not return targets, pull them from the literature
    if drug_genes is None:

        # Search EBC for a drug target, via binding annotations
        drug_gene_list = EBC_api.query_drug_target(drug)

    else:
        # If targets are from Pharos, map them to their Uniprot IDs
        drug_gene_list = []
        for gene in drug_genes:
            drug_gene_list.append(EBC_api.resolve_gene_to_EntrezGeneID(gene))

    # If either disease or drug list comes up empty
    # the query has failed and we return a statement to that effect

    if dis_gene_list is None and drug_gene_list is None:
        return "Drug and disease not recognized, better luck next time..."
    elif dis_gene_list is None:
        return "Disease not recognized"
    elif drug_gene_list is None:
        return "Drug not recognized"

    # If we have targets
    else:

        PMIDs = EBC_api.query_chemical_disease(drug, disease, get_PMIDs=True)
        print(PMIDs)

        # Select the top 25 genes from the disease gene list for GO enrichment
        dis_genes = [[EBC_api.resolve_EntrezGeneID_to_NCBIGeneName(x),x] for x in dis_gene_list]

        dis_genes = pd.DataFrame(dis_genes, columns=["Gene", "Entrez ID"])
        dis_gene_list = list(map(int, dis_gene_list))

        # dis_gene_names = EBC_api.get_disease_gene_list(disease, freq_correct=True, gene_names=True)
        # print(dis_gene_names[:10])
        # disease_tissues = pharos_api.get_tissues_oi(dis_gene_names[:5])
        # print(disease_tissues)
        # print(drug_genes)

        # Get drug tissues
        # drug_tissues = pharos_api.get_tissues_oi(drug_genes)
        # print(drug_tissues)



        drug_genes = [[EBC_api.resolve_EntrezGeneID_to_NCBIGeneName(x), x] for x in drug_gene_list]
        drug_genes = pd.DataFrame(drug_genes, columns=["Gene", "Entrez ID"])
        
        # Get the GO terms for the drug targets
        drug_gene_list = list(map(int,drug_gene_list))

        drug_targets = GO_API.get_GO_terms(drug_gene_list)

#        Get GO Enrichment statistics
        go_result = GO_API.calc_GO_enrichment(dis_gene_list, os.path.join(results_dir, out_name), target_list=drug_targets)
        if output_full is False:
            go_result = go_result.loc[go_result['rejected'] == 1.0, ['namespace', 'name', 'p', 'q']]
            go_result = go_result.sort(['q'])

        result = {"GOENRICH":go_result, "drug_genes":drug_genes, "disease_genes":dis_genes}

        if gen_interaction_image:
            subprocess.check_call(['dot', '-Tpng', os.path.join(results_dir, out_name)+ '.dot', '-o', os.path.join(results_dir, out_name) + '.png'])
            result["image_file"] = os.path.join(results_dir, out_name) + '.png'
        if gen_tissues_image:
            pass



    # return(drug_tissues)
    return(result)

# unit testing
if __name__ == "__main__":
    drug="lisinopril"

    disease="hypertension"

    GO_API = go_api.GO_api("./ncats/GO_api/GO_DB")

    print(drug,disease)
    result = Q2_query(drug, disease)
    # if type(result) is not str:
    #     print(result)
    # else:
    #     print(result)

    drug="tacrine"

    disease="alzheimer-disease"

    print(drug,disease)
    result = Q2_query(drug, disease)
    # if type(result) is not str:
    #     print(result)
    # else:
    #     print(result)
