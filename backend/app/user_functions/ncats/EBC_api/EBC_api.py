import pickle
import numpy as np
from fuzzywuzzy import process
import os, sys

# Path to the sys
EBC_path = os.path.abspath(os.path.dirname(__file__))

EBC_reference_dir = os.path.join(EBC_path, "EBC_DB")

# Import reference dictionary objects
pkl_file = open(os.path.join(EBC_reference_dir, 'chem_GNDT_dict.pkl'), 'rb')
chem_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir, 'disease_GNDT_dict.pkl'), 'rb')
disease_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'disease_matching_dict.pkl'), 'rb')
disease_matching_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'gene_dict.pkl'), 'rb')
gene_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'gene_idf.pkl'), 'rb')
gene_idf = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'chemical_disease_GNDT_dict.pkl'), 'rb')
chem_disease_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'chemical_gene_GNDT_dict.pkl'), 'rb')
chem_gene_dict = pickle.load(pkl_file)
pkl_file.close()
#
pkl_file = open(os.path.join(EBC_reference_dir,'disease_gene_GNDT_dict.pkl'), 'rb')
dis_gene_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(EBC_reference_dir,'gene_gene_GNDT_dict.pkl'), 'rb')
gene_gene_dict = pickle.load(pkl_file)
pkl_file.close()

"""
query_drug_target

Querys a drug for a given relationship with a gene, default is "B", a binding relationships

Input:
drug - string of the drug common name
relationship - desired disease relationship
gene_names - return a list of gene names instead of ENTREZ gene IDs
n - max number of genes to return

Output:
ordered_genes - returns a list of the top n genes ordered by the frequency of that genes annotation to the drug
If no gene targets are known or the drug isn't in our database, it returns None
"""
def query_drug_target(drug, relationship="B", gene_names=False, n=25):
        chem_id = chem_dict.get(drug.lower())

        count = dict()
        if chem_id in chem_gene_dict:
            for entry in chem_gene_dict.get(chem_id):
                id = entry[0]
                if id in gene_dict and id in gene_idf:
                    if entry[1] == relationship:
                        if id in count:
                            count[id] += 1
                        else:
                            count[id] = 1

            freq_corrected_genes = dict()
            if gene_names:
                for gene in count:
                    freq_corrected_genes[gene_dict.get(gene)] = np.multiply(float(count.get(gene)), gene_idf.get(gene))
            else:
                for gene in count:
                    freq_corrected_genes[gene] = np.multiply(float(count.get(gene)), gene_idf.get(gene))
            ordered_genes = sorted(freq_corrected_genes, key=lambda key: freq_corrected_genes[key], reverse=True)
            if len(ordered_genes) > 0:
                return(ordered_genes[:n])
            else:
                return(None)

        else:
            return(None)

"""
query_drug_target

Querys a drug for a given relationship with a gene, default is "B", a binding relationships

Input:
drug - string of the drug common name
relationship - desired disease relationship
gene_names - return a list of gene names instead of ENTREZ gene IDs
n - max number of genes to return

Output:
ordered_genes - returns a list of the top n genes ordered by the frequency of that genes annotation to the drug
If no gene targets are known or the drug isn't in our database, it returns None
"""
def query_disease(disease, relationship="U", gen_counts=False):
    disease_id = disease_dict.get(disease.lower())
    disease_set = set()
    count = dict()
    if disease_id in dis_gene_dict:
        for entry in dis_gene_dict.get(disease_id):
            id = entry[0]
            if id in gene_dict:
                if entry[1] == relationship:
                    if gen_counts:
                        if id in count:
                            count[id] += 1
                        else:
                            count[id] = 1
                    else:
                        disease_set.add(id)
        if gen_counts:
            return(count)
        else:
            return(list(disease_set))
    else:
        return(None)

################################

def query_chemical_disease(drug, disease, relationship="T", gen_counts=False, get_PMIDs=False):
    chem_id = chem_dict.get(drug.lower())
    disease_id = disease_dict.get(disease.lower())
    disease_set = set()
    if gen_counts:
        count = {'Mp': 0, 'C': 0, 'T': 0, 'J': 0, 'Sa': 0, 'Pr': 0, 'Pa': 0}
        if chem_id in chem_disease_dict:
            for entry in chem_disease_dict.get(chem_id):
                id = entry[0]
                if id == disease_id:
                    # print(drug, disease, entry[1])
                    if gen_counts:
                        count[entry[1]] += 1
                    # if entry[1] == relationship:
                    #     if gen_counts:
                    #         if id in count:
                    #             count[id] += 1
                    #         else:
                    #             count[id] = 1
                    #     else:
                    #         disease_set.add(id)
        return(count)

    elif get_PMIDs:
        PMIDs = []
        if chem_id in chem_disease_dict:
            for entry in chem_disease_dict.get(chem_id):
                id = entry[0]
                # prin t(entry)
                if id == disease_id:
                    if entry[1] == relationship:
                        PMIDs.append(entry[2])
        return(PMIDs)
        # if gen_counts:
        #     return(count)
        # else:
        #     return(list(disease_set))
    else:
        return(None)

def query_gene(gene, relationship="B", gene_name=False, gen_counts=False):
    gene_set = set()
    count = dict()
    for entry in gene_gene_dict.get(gene):
        if entry[1] in relationship:
            if gene_name:
                gene_set.add(gene_dict.get(entry[0]))
            else:
                gene_set.add(entry[0])
    return(list(gene_set))

def resolve_gene_to_EntrezGeneID(geneName):
    return(gene_dict.get(geneName))

def resolve_EntrezGeneID_to_NCBIGeneName(geneName):
    return(gene_dict.get(geneName))

def get_disease_gene_list(disease, relationship = "U", gene_names = False, freq_correct = False, return_scores=False):

    disease_id = disease_dict.get(disease.lower())
    if disease_id in dis_gene_dict:
        count = dict()
        for entry in dis_gene_dict.get(disease_id):
            id = entry[0]
            if id in gene_dict and id in gene_idf:
                if entry[1] == relationship:
                    if id in count:
                        count[id] += 1
                    else:
                        count[id] = 1
        if freq_correct:
            freq_corrected_genes = dict()
            if gene_names:
                for gene in count:
                    freq_corrected_genes[gene_dict.get(gene)] = np.multiply(float(count.get(gene)), gene_idf.get(gene))
            else:
                for gene in count:
                    freq_corrected_genes[gene] = np.multiply(float(count.get(gene)), gene_idf.get(gene))
            if return_scores:
                return(freq_corrected_genes)
            else:
                order_genes2 = sorted(freq_corrected_genes, key=lambda key: freq_corrected_genes[key], reverse=True)

            return(order_genes2)
        else:
            order_genes = sorted(count, key=lambda key: count[key], reverse=True)
            return(order_genes)
    else:
        return(None)


# Can be generalized to any desired relationship
def query_for_drug_disease_relationships(filePath, outPath, themes= ['T', 'C', 'Mp', 'J', 'Sa', 'Pr', 'Pa']):
    drugs = []
    with open(filePath) as inFile:
        for line in inFile.readlines():
            drugs.append(line.strip())

    with open(outPath, "w+") as outFile:
        outFile.write("drug,disease," + ",".join(themes) + "\n")
        for drug in drugs:
            for disease in diseases:
                outFile.write(",".join([drug, disease]))
                print(drug, disease)
                results = query_chemical_disease(drug, disease, gen_counts=True)
                if results is None:
                    outFile.write(",-,-,-,-,-,-")
                else:
                    for theme in themes:
                        outFile.write("," + str(results.get(theme)))
                outFile.write("\n")

def query_term_for_matches(queryTerm, searchSet = disease_matching_dict, numMatches = 5):

    if "," in queryTerm:
        search1_terms = queryTerm.split(", ")
        quick_search_term = " ".join(search1_terms[::-1])
        if quick_search_term in disease_dict:
            return(quick_search_term, 0)
    else:
        search1_terms = queryTerm.split(" ")

    # Use partial keys to narrow search space
    newSearchSpace = set()
    for partial_term in search1_terms:
        if partial_term in searchSet:
            newSearchSpace = newSearchSpace.union(searchSet.get(partial_term))
    if len(newSearchSpace) == 0:
        return(None, -1)

    results = process.extract(queryTerm, newSearchSpace, limit=numMatches)
    return(results, 1)


if __name__ == "__main__":
    pkl_file = open('./chem_GNDT_dict.pkl', 'rb')
    chem_dict = pickle.load(pkl_file)
    pkl_file.close()

    pkl_file = open('./disease_GNDT_dict.pkl', 'rb')
    disease_dict = pickle.load(pkl_file)
    pkl_file.close()

    pkl_file = open('./disease_matching_dict.pkl', 'rb')
    disease_matching_dict = pickle.load(pkl_file)
    pkl_file.close()

    pkl_file = open('./gene_GNDT_dict.pkl', 'rb')
    gene_dict = pickle.load(pkl_file)
    pkl_file.close()
    with open("../disease_failues_matched.txt", "w+") as outFile:
        with open("../disease_failures.txt", "r") as inFile:
            for disease in inFile.readlines():
                disease = disease.strip().lower()
                print(disease)
                matches, method = query_term_for_matches(disease, disease_matching_dict)
                if method == 0:
                    match = matches
                    new_key = disease_dict.get(matches)
                elif method == -1:
                    new_key = "-"
                    match = "-"
                else:
                    match = matches[0][0]
                    new_key = disease_dict.get(matches[0][0])

                print(disease + "\t" + match + "\t" + new_key + "\t" + str(method) + "\n")
                outFile.write(disease + "\t" + match + "\t" + new_key + "\t" + str(method) + "\n")

if __name__ == "__main__":

    # INSERT UNIT TEST HERE
    exit()

