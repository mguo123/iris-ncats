import pickle
import numpy as np
import os


########## LOAD REFERENCE FILES INTO MEMORY FOR FASTER QUERIES ##########

# to Ncats path
GNBR_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

GNBR_reference_dir = os.path.join(GNBR_path, "DB_data/GNBR_DB")

# Python dictionaries containing different term mappings, resolve common language keys to MEDIC Identities.
pkl_file = open(os.path.join(GNBR_reference_dir, 'chem_GNBR_dict.pkl'), 'rb')
chem_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir, 'CTD_chem_mapping_dict.pkl'), 'rb')
chem_entity_resolution = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir, 'disease_GNBR_dict.pkl'), 'rb')
disease_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir,'gene_dict.pkl'), 'rb')
gene_dict = pickle.load(pkl_file)
pkl_file.close()

# Dictionary of different substrings for partial matching and searching
pkl_file = open(os.path.join(GNBR_reference_dir,'disease_matching_dict.pkl'), 'rb')
disease_matching_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir,'gene_idf.pkl'), 'rb')
gene_idf = pickle.load(pkl_file)
pkl_file.close()

# Dictionary objects of MEDIC identities and the GNBR theme relationships between them
pkl_file = open(os.path.join(GNBR_reference_dir,'chemical_disease_GNBR_dict.pkl'), 'rb')
chem_disease_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir,'chemical_gene_GNBR_dict.pkl'), 'rb')
chem_gene_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir,'disease_gene_GNBR_dict.pkl'), 'rb')
dis_gene_dict = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open(os.path.join(GNBR_reference_dir,'gene_gene_GNBR_dict.pkl'), 'rb')
gene_gene_dict = pickle.load(pkl_file)
pkl_file.close()


def query_drug_target(drug, relationship="B", gene_names=False, n=25):
    """
    query_drug_target

    Querys a drug for a given relationship with a gene, default is "B", a binding relationships

    :param drug: string of the drug common name
    :param relationship: desired drug-gene relationship
    :param gene_names: return a list of gene names instead of ENTREZ gene IDs
    :param n: max number of genes to return
    :return: returns a list of the top n genes ordered by the frequency of that genes annotation to the drug
    If no gene targets are known or the drug isn't in our database, it returns None
    """
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


def query_disease(disease, relationship="U", gen_counts=False):
    """
    query_disease

    Querys a disease for a given relationship with a gene, default is "U", a causal mutation relationship


    :param disease: string of the common disease name
    :param relationship: desired disease relationship
    :param gen_counts: Whether a dictionary with the counts of each relationship instance should be return or not
    :return: disease_set - returns a list or dictionary of the genes annotated to the disease with or without counts
    If no gene interactions are known or the disease isn't in our database, it returns None
    """
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
            if len(disease_set.keys()) == 0:
                return (None)
            else:
                return(list(disease_set))
    else:
        return(None)

def query_chemical_disease(drug, disease, relationship="T", gen_counts=False, get_PMIDs=False):
    """
    :param drug: string of the common drug name
    :param disease: string of the common disease name
    :param relationship: desired disease relationship
    :param gen_counts: return a dictionary object containing a count of each themes occurrence for the drug and disease
    :param get_PMIDs: list of PMIDs support the indicated relationship
    :return:
    """
    chem_id = chem_dict.get(drug.lower())
    disease_id = disease_dict.get(disease.lower())
    disease_set = set()
    if gen_counts:
        count = {'Mp': 0, 'C': 0, 'T': 0, 'J': 0, 'Sa': 0, 'Pr': 0, 'Pa': 0}
        if chem_id in chem_disease_dict:
            for entry in chem_disease_dict.get(chem_id):
                id = entry[0]
                if id == disease_id:
                    if gen_counts:
                        count[entry[1]] += 1
        return(count)

    elif get_PMIDs:
        PMIDs = []
        if chem_id in chem_disease_dict:
            for entry in chem_disease_dict.get(chem_id):
                id = entry[0]
                if id == disease_id:
                    if entry[1] == relationship:
                        PMIDs.append(entry[2])
        return(PMIDs)
    else:
        return(None)

def query_gene(gene, relationship="B", gene_name=False, gen_counts=False):
    """
    :param gene: the NCBI gene name
    :param relationship: desired disease relationship
    :param gene_name: return a list of NCBI gene names instead of ENTREZ gene IDs
    :return: gene_set - returns a list of with all genes that have the indicated relationship
    """
    gene_set = set()
    for entry in gene_gene_dict.get(gene):
        if entry[1] in relationship:
            if gene_name:
                gene_set.add(gene_dict.get(entry[0]))
            else:
                gene_set.add(entry[0])
    return(list(gene_set))

def resolve_gene_to_EntrezGeneID(geneName):
    """
    resolve_gene_to_EntrezGeneID

    Querys an NCBI gene name and returns the Entrez Gene ID

    :param geneName: NCBI gene name
    :return: Entrez Gene ID (str)
    """
    return(gene_dict.get(geneName))

def resolve_EntrezGeneID_to_NCBIGeneName(geneName):
    """
    resolve_EntrezGeneID_to_NCBIGeneName

    Querys an Entrez Gene ID and returns the NCBI gene name

    :param geneName: Entrez gene ID (str)
    :return: NCBI gene name (str)
    """
    return(gene_dict.get(geneName))

def get_disease_gene_list(disease, relationship = "U", gene_names = False, freq_correct = False, return_scores=False):
    """
    :param disease: common disease name
    :param relationship: GNBR relationship
    :param gene_names: return NCBI gene names instead of EntrezIDs
    :param freq_correct: Correct for the frequency of annotation
    :param return_scores: Return the score associated with the frequency correction
    :return: list or dictionary of genes with the desired relationship to the disease
    """
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
        if len(count.keys()) == 0:
            return None

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

def query_term_for_matches(queryTerm, searchSet = disease_matching_dict, numMatches = 5):
    """
    query_term_for_matches

    Queries a common name for a disease to identify potential matches that are in the database.
    Essentially a "Did you mean?" kind of functionality.

    :param queryTerm: common name for partial matching search
    :param searchSet: Python dictionary containg partial strings and associated search spaces
    :param numMatches: number of matches to return, sorted by similarity score
    :return: results, lise of potential matches
    """
    from fuzzywuzzy import process

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


def get_MEDICID(term):
    """
    :param term: common name for the term
    :return: MEDIC ID for that term
    """
    if term in disease_dict:
        return disease_dict.get(term)
    elif term in chem_dict:
        return chem_dict.get(term)
    else:
        return 'ID:Unknown'

def query_treatments_for_disease(disease, relationship="T"):
        """
        :param disease: common disease name
        :param relationship: GNBR relationship
        :param gene_names: return NCBI gene names instead of EntrezIDs
        :param freq_correct: Correct for the frequency of annotation
        :param return_scores: Return the score associated with the frequency correction
        :return: list or dictionary of genes with the desired relationship to the disease
        """
        disease_id = disease_dict.get(disease.lower())

        treatments = dict()
        if disease_id in chem_disease_dict:
            count = dict()
            for entry in chem_disease_dict.get(disease_id):
                print(entry)
                id = entry[0]
                if id in chem_dict:
                    print("hit")
                    if entry[1] == relationship:
                        if id in count:
                            count[id] += 1
                        else:
                            count[id] = 1

            if len(count.keys()) == 0:
                return None
        count_sorted = sorted(count, key=lambda key: count[key], reverse=True)
        out_list = []
        for key in count_sorted:
            out_list.append([chem_entity_resolution.get(key),key, count.get(key)])
        print(out_list)

if __name__ == "__main__":

    # print(query_term_for_matches("sickle cell trait"))
    print(query_treatments_for_disease("rheumatoid arthritis"))








