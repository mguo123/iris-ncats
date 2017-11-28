import os
import pickle
import sys
import subprocess

from optparse import OptionParser

import pandas as pd
from pubmed_lookup import PubMedLookup, Publication

from app.user_functions.ncats.Q2.TiGER_api import TiGER_api
from app.user_functions.ncats.Q2.GNBR_api import GNBR_api
from app.user_functions.ncats.Q2.Pharos_api import pharos_api
from app.user_functions.ncats.Q2.GO_api import go_api

# from GNBR_api import GNBR_api
# from Pharos_api import pharos_api
# from GO_api import go_api
# from TiGER_api import TiGER_api


############################################# DEFINE PATHS ##############################################################


ncats_path = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))


GO_API = go_api.GO_api(os.path.join(ncats_path, "DB_data/GO_DB"))


############################################# Functions ##############################################################

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
def get_PMID(PMID):
    email = ''
    url = 'http://www.ncbi.nlm.nih.gov/pubmed/' + PMID

    try:

        lookup = PubMedLookup(url, email)
        publication = Publication(lookup)
        return(publication.title)
    except:
        return "Title not found"

def Q2_query(QDrug, QDisease, options):

    # Pre-process query
    drug = QDrug.strip().lower()

    disease = QDisease.strip().lower()


    # get MEDIC ids
    drug_id = GNBR_api.get_MEDICID(drug)
    disease_id = GNBR_api.get_MEDICID(disease)


    # Generate prefix for output file, is drug_disease
    out_name = drug + "_" + disease.replace(" ", "-").lower()

    if options.verbose and options.batchFile is None:
        print("Querying GNBR for disease genes...")
    # Get list of genes causally annotated to a disease
    dis_gene_list = GNBR_api.get_disease_gene_list(disease, freq_correct=True)

    # If first query did not work, try getting nearest matches for the query and try again
    if dis_gene_list is None:

        # Get matches
        disease2, match_type = GNBR_api.query_term_for_matches(disease)

        # If not fuzzy matching, but simple order matching, use the new disease query and proceed
        if match_type == 0:
            dis_gene_list = GNBR_api.get_disease_gene_list(disease2, freq_correct=True)

    if options.verbose and options.batchFile is None:
        print("Querying Pharos for drug targets...")
    # Get list of drug targets from Pharos
    drug_genes = pharos_api.get_ligand_targets(drug)

    # If Pharos did not return targets, pull them from the literature
    if drug_genes is None:
        if options.verbose and options.batchFile is None:
            print("Pharos did not contain drug target information, querying GNBR for drug targets...")
        # Search GNBR for a drug target, via binding annotations
        drug_gene_list = GNBR_api.query_drug_target(drug)

    else:
        # If targets are from Pharos, map them to their Uniprot IDs
        drug_gene_list = []
        print("If targets are from Pharos, map them to their Uniprot IDs")
        for gene in drug_genes:
            drug_gene_list.append(GNBR_api.resolve_gene_to_EntrezGeneID(gene))

    # If either disease or drug list comes up empty
    # the query has failed and we return a statement to that effect

    if dis_gene_list is None and drug_gene_list is None:
        print("ERROR: Drug and disease not recognized")
        return "ERROR: Drug and disease not recognized"

    elif dis_gene_list is None:
        print("ERROR: Disease not recognized")
        return "ERROR: Disease not recognized"

    elif drug_gene_list is None:
        print("ERROR: No drug targets found")
        return "ERROR: No drug targets found"

    # If we have targets
    else:
        print("drug (%s) and disease (%s) found. Processing..." % (drug, disease) )
        # Generate output directory
        overall_path = os.path.abspath(os.path.dirname(__file__))
        results_dir = os.path.join(*[overall_path, options.outPath, out_name])

        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        # Select the top 25 genes from the disease gene list for GO enrichment

        print("Getting disease genes Calling resolve_EntrezGeneID_to_NCBIGeneName")
        dis_genes = pd.DataFrame([[GNBR_api.resolve_EntrezGeneID_to_NCBIGeneName(str(x)),x] for x in dis_gene_list], columns=["Gene", "Entrez ID"])
        dis_genes_short = dis_genes[:min(len(dis_genes), 5)]
        dis_gene_list = list(map(int, dis_gene_list))

 
        print("Getting drug genes Calling resolve_EntrezGeneID_to_NCBIGeneName")
        drug_genes = [[GNBR_api.resolve_EntrezGeneID_to_NCBIGeneName(x), x] for x in drug_gene_list]
        drug_genes = pd.DataFrame(drug_genes, columns=["Gene", "Entrez ID"])
        drug_genes_short = drug_genes[:min(5, len(drug_genes))]

        # # Get tissue information
        # print("Getting tissue information resolved to drug via: resolve_EntrezGeneID_to_NCBIGeneName")
        # tissue_df_drug = TiGER_api.get_tissue_counts([GNBR_api.resolve_EntrezGeneID_to_NCBIGeneName(str(x)) for x in drug_gene_list])
        # tissue_df_drug_short = tissue_df_drug[:min(5, len(tissue_df_drug))]

        # Get the GO terms for the drug targets
        drug_gene_list = list(map(int,drug_gene_list))

        drug_targets = GO_API.get_GO_terms(drug_gene_list)

#        Get GO Enrichment statistics
        print("Getting Go Enrichment statistics")
        if options.gen_image:
            go_result = GO_API.calc_GO_enrichment(dis_gene_list, os.path.join(results_dir, out_name), target_list=drug_targets, gen_image=True)
        else:
            go_result = GO_API.calc_GO_enrichment(dis_gene_list, "",target_list=drug_targets)

        go_result['gene_target'] = go_result['term'].isin(drug_targets)

        go_result = go_result.loc[go_result['rejected'] == 1.0, ['name', 'term', 'p', 'q', 'gene_target']]
        go_result = go_result.sort_values(by=['gene_target', 'q'], ascending=[False, True])
        go_result.to_csv(os.path.join(results_dir, out_name + "_GO_pathway_enrichment.csv"), mode="w+", index_label=False, index=False)


        # Get GO Enrichment statistics
        go_result_short = go_result[:min(5, len(go_result))]

        # Start saving results
        result = {"GOENRICH":go_result, "drug_genes":drug_genes, "disease_genes":dis_genes, "drug_id": drug_id, "disease_id":disease_id,
                  "GOENRICH_short":go_result_short, "drug_genes_short":drug_genes_short, "disease_genes_short":dis_genes_short,
                  }


        # Get tissue information
        print("Getting tissue information resolved to disease via: resolve_EntrezGeneID_to_NCBIGeneName")
        tissue_df_dis = TiGER_api.get_tissue_counts([GNBR_api.resolve_EntrezGeneID_to_NCBIGeneName(str(x)) for x in dis_gene_list])
        if tissue_df_dis is not None:
            tissue_df_dis_short = tissue_df_dis[:min(5, len(tissue_df_dis))]
            result["tissue_df_dis"] = tissue_df_dis
            result["tissue_df_dis_short"] = tissue_df_dis_short

        print('Generating Image')
        if options.gen_image:
            file_name = os.path.join(results_dir, out_name + '.png')
            if os.path.exists(os.path.join(results_dir, out_name + '.dot')):
                subprocess.check_call(['dot', '-Tpng', os.path.join(results_dir, out_name + '.dot'), '-o', file_name])

            result["image_file"] = file_name

        # Get Pubmed id
        print("Getting Pubmed IDs")
        if options.gen_pubmed:
            PMIDs = GNBR_api.query_chemical_disease(drug, disease, get_PMIDs=True)
            len(PMIDs)
            # print(PMIDs.shape)
            # print('Saving pubmed PMIDs to ', results_dir, out_name,"_PMIDs.csv")
            # PMID_df.to_csv(os.path.join(results_dir, out_name + "_PMIDs.csv"), mode="w+",
            #                      index_label=False, index=False, header=False)
            print("Getting Pubmed Titles")
            if len(PMIDs) > 0:
                ################### THIS IS JUST TOP 10 FOR NOW, SEEMS LIKE THE API CALLS TAKE SOME TIME
                PMID_df = pd.DataFrame([[x, get_PMID(x)] for x in PMIDs[:min(10, len(PMIDs))]], columns=["PMIDS", "Title"])
                
                # will show top 5
                PMID_df_short = PMID_df[:min(5, len(PMID_df))]
                result["pubmed"] = PMID_df
            
                result["pubmed_short"] = PMID_df_short

                print('Saving pubmed PMIDs to ', results_dir, out_name,"_PMIDs.csv")
                PMID_df.to_csv(os.path.join(results_dir, out_name + "_PMIDs.csv"), mode="w+",
                                 index_label=False, index=False, header=False)
            
            else:
                result["pubmed"] = 'no PMIDs found'

    # return(drug_tissues)
    return(result)

# unit testing
if __name__ == '__main__':
    usage = "usage: \n%prog [options] <common drug name> <common_disease_name>\n\nOR\n\n%prog -f <path_to_batch_file>"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", dest="batchFile", type="str", help="path to batch query file")
    parser.add_option("-v", action="store_true", help="verbose", dest="verbose", default=False)
    parser.add_option("-g", action="store_true", help="generate graphic", dest="gen_image", default=False)
    parser.add_option("-p", action="store_true", help="output pubmed IDs", dest="gen_pubmed", default=False)
    parser.add_option("-o", dest="outPath", help="path to output directory", default="results/")

    (options, args) = parser.parse_args()

    if options.batchFile is None:
        if len(args) < 2:
            print("ERROR: NO DRUG AND/OR DISEASE GIVEN\n")
            print(parser.print_help())
            exit()
        else:
            print("Initializing...")
            overall_path = os.path.abspath(os.path.dirname(__file__))
            GO_API = go_api.GO_api(os.path.join(overall_path, "../DB_data/GO_DB"))
            drug = args[0].lower()
            disease = " ".join(args[1:]).lower()
            if options.verbose:
                print("Query = " + drug + " and " + disease)
            print("Running Query...")
            Q2_query(drug, disease, options)
    else:
        print("Running Batch Query")
        overall_path = os.path.abspath(os.path.dirname(__file__))
        GO_API = go_api.GO_api(os.path.join(overall_path, "../DB_data/GO_DB"))
        with open(options.batchFile, "r") as inBatch:
            for line in inBatch:
                info = line.strip().split("\t")
                if len(info) != 2:
                    print("Poorly formed query...")
                    print("ERROR IN LINE: " + line)
                    exit()

                drug = info[0].lower()
                disease = info[1].lower()
                if options.verbose:
                    print("Query= " + " and ".join([drug, disease]))
                Q2_query(drug, disease, options)


