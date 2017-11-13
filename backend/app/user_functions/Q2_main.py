from ncats.EBC_api import EBC_api
from ncats.Pharos_api import pharos_api
from ncats.GO_api import go_api

import pickle
import subprocess
import os

results_dir = "./ncats/results/Q2/"

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

def Q2_query(QDrug, QDisease, gen_image=False, output_full=False):

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
        # Select the top 25 genes from the disease gene list for GO enrichment
        dis_gene_list = list(map(int, dis_gene_list))[:25]

        # Get the GO terms for the drug targets
        drug_gene_list = list(map(int,drug_gene_list))
        drug_targets = go_api.get_GO_terms(drug_gene_list)

        # Get GO Enrichment statistics
        result = go_api.calc_GO_enrichment(dis_gene_list, os.path.join(results_dir, out_name), target_list=drug_targets)
        if gen_image:
            subprocess.check_call(['dot', '-Tpng', os.path.join(results_dir, out_name)+ '.dot', '-o', os.path.join(results_dir, out_name) + '.png'])

    return(result)

if __name__ == "__main__":
    drug="lisinopril"

    disease="hypertension"

    print(drug,disease)
    result = Q2_query(drug, disease)
    if type(result) is not str:
        print(result.shape)
    else:
        print(result)
