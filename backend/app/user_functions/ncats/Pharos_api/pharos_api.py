import requests
import nltk
import os
## REPLACE WITH YOUR LOCAL PATH
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # should yield the same as .. (with subfolders rcsc, scripts, and results)
RSCS_DIR = os.path.join(PARENT_DIR, 'shared_data')
PATH_TO_GIANT_TISSUES_LIST=os.path.join(RSCS_DIR, 'list_of_GIANT_tissues.txt')


def query(data_type, uid):
    """
    query

    Generates a Pharos API query and returns the results

    Input:
    data_type - the Pharos data type (e.g. ligand, target ..)
    uid - Pharos concept identifier
    Output:
    JSON containing the query results
    """
    url = "https://pharos.nih.gov/idg/api/v1/{}({})?view=full".format(data_type, uid)
    try:
        resp = requests.get(url).json()
    except:
        resp = None
    return resp


####################################################################################################
####################################################################################################
###################################    DISEASE TARGET FUNCTIONS   ##################################
####################################################################################################
####################################################################################################
def get_disease_id(disease_common_name):
    """
    get_disease_id

    Queries Pharon API and return the Pharon ID for the disease if it is in Pharos

    Input: 
    disease_common_name - name of the disease query, will search for this term

    output:
    disease_id = Pharoid uid for the queried disease, if none exists, return None



    """
    disease_common_name = disease_common_name.replace(" ", "%20").lower()
    p = nltk.PorterStemmer()
    disease_common_name = p.stem(p.stem(disease_common_name))
    url = "https://pharos.nih.gov/idg/api/v1/diseases/search?q={}&top=10".format(disease_common_name)
    resp = requests.get(url).json()
    if len(resp["content"]) != 0:
        disease_id = resp["content"][0]["id"]
    else:
        disease_id = None

    return disease_id

def disease_synonym_data(qID):
    url = "https://pharos.nih.gov/idg/api/v1/diseases/{}/synonyms".format(qID)
    info = requests.get(url).json()
    return info

####################################################################################################
####################################################################################################
###################################        DRUG FUNCTIONS         ##################################
####################################################################################################
####################################################################################################

def get_ligand_id(ligand_common_name):
    """
    get_ligand_id

    Queries the Pharos API and returns the Pharos ID for the ligand if it is in Pharos

    Input:
    ligand_common_name - the common name for the drug, this will be the search term for Pharos

    Output:
    ligand_id - Pharos uid for the queried drug, if one exists, else None is returned
    """
    url = "https://pharos.nih.gov/idg/api/v1/ligands/search?q={}&top=10".format(ligand_common_name)

    resp = requests.get(url).json()
    if len(resp["content"]) != 0:
        ligand_id = resp["content"][0]["id"]
    else:
        ligand_id = None
    return ligand_id


def get_ligand_targets(ligand):
    """
    get_ligand_targets

    Queries the Pharos API and returns the drug targets for the ligand provided
    Input:
    ligand - Pharos uid for the drug

    Output:
    putative_targets - a list of putative targets, return None if there were no known targets in Pharos
    """
    qID = get_ligand_id(ligand)
    resp = query("ligands", qID)
    if resp is None:
        return resp
    else:
        putative_targets = []
        for link in resp["links"]:
            if link["kind"] == "ix.idg.models.Target":
                for prop in link["properties"]:
                    if prop["label"] == 'IDG Target':
                        putative_targets.append(prop['term'])

    return putative_targets


####################################################################################################
####################################################################################################
###################################     GENE TARGET FUNCTIONS     ##################################
####################################################################################################
####################################################################################################
def get_target_id(target_common_name):
    """
    get_target_id(target_common_name)

    Input: 
    target_common_name - drug name

    Output: 
    target_id - corresponding Pharos ID for the target
    """
    try:
        resp = query("targets", target_common_name)
        target_id = resp['id']
    except:
        target_id = None

    return target_id



def target_prop_data(qID, label):

    """
    target_prop_data(qID, label)

    Input: 
    qID - Pharos ID for the gene target
    label - property of interest

    Output: 
    info - list of dictionaries containing the properties data with the specified label
    """
    url = "https://pharos.nih.gov/idg/api/v1/targets/{}/properties(label={})".format(qID, label)
    info = requests.get(url).json()
    return info

def target_link_data(qID, kind):

    """
    target_link_data(qID, label)

    Input: 
    qID - Pharos ID for the gene target
    label - link of interest

    Output: 
    info - list of dictionaries containing the link data with the specified label
    """
    url = "https://pharos.nih.gov/idg/api/v1/targets/{}/links(kind={})".format(qID, kind)
    info = requests.get(url).json()
    return info

def target_synonym_data(qID):
    url = "https://pharos.nih.gov/idg/api/v1/targets/{}/synonyms".format(qID)
    info = requests.get(url).json()
    return info

def get_target_interactions(qID):
    for synonym in target_synonym_data(qID):
        if synonym["label"] == "UniProt Accession":
            uniID = synonym["href"][::-1].split("/")[0][::-1]
            url = "https://www.ebi.ac.uk/proteins/api/proteins/interaction/{}".format(uniID)
            return(requests.get(url).json())


def get_pathways(target_name):
    """
    get_pathways(target_name)

    Input: 
    targets - drug name

    Output: 
    pathways, list of pathways relevant for disease (use KEGG Pathway, Reactome Pathway, PathwayCommons*, WikiPathways Pathway)
    """
    qID = get_target_id(target_name)

    # Get tissue information from the properties for that target
    pathway_info = target_prop_data(qID, "*Pathway")

    # Grab PATHWAY annotations for each gene
    pathways = set()
    for prop in pathway_info:
        pathways.add(prop['term'])
    if len(pathways) == 0:
        pathways.add("unknown")
    
    return pathways


def get_GO_terms(target_name):
    """
    get_GO_terms(target_name)

    Input: 
    targets - drug name

    Output: 
    GO_terms - list of GO process, function, and component (list of 3 dictionaries) terms for that drug
    """
    # Query for the Pharos ID of the Gene target
    qID = get_target_id(target_name)

    # Get tissue information from the properties for that target
    GO_info = target_prop_data(qID, "*GO")

    # Grab GO term annotations for each gene
    GO_process = set()
    GO_function = set()
    GO_component = set()
    for prop in GO_info:
        if prop["label"] == 'GO Process':
            GO_process.add(prop['term'])
        elif prop["label"] == 'GO Function':
            GO_function.add(prop['term'])
        elif prop["label"] == 'GO Component':
            GO_component.add(prop['term'])
    if len(GO_process) == 0:
        GO_process.add("unknown")
    if len(GO_function) == 0:
        GO_function.add("unknown")
    if len(GO_component) == 0:
        GO_component.add("unknown")

    return GO_process, GO_function, GO_component


def get_target_tissue(target_name):
    """
    get_target_tissue(target_name)

    Input: 
    targets - drug name

    Output: 
    tissues - list of tissues for that drug
    """
    # Query for the Pharos ID of the Gene target
    qID = get_target_id(target_name)
    if qID is None:
        return []
    # Get tissue information from the properties for that target
    tissue_info = target_prop_data(qID, "*Tissue")

    # Grab UniProt Tissue annotations for each gene
    tissues = set()
    for prop in tissue_info:
        if prop["label"] == 'UniProt Tissue':
            tissues.add(prop['term'])
    if len(tissues) == 0:
        tissues.add("unknown")

    return tissues


def get_tissues_oi(targets):
    """
    get_tissues_oi(targets)

    Input: 
    targets - list of gene targets

    Output: 
    tissues_oi - list of tissues where those genes are expressed
    """
    tissues_oi = dict()
    # Get tissues associated with each target and reformat for compatibility with GIANT
    for target in targets:
        for x in get_target_tissue(target):
            tissue = x.lower().replace(" ", "_")
            if tissue in tissues_oi:
                tissues_oi[tissue] +=1
            else:
                tissues_oi[tissue] = 1
    return tissues_oi

"""
init_GIANT_queries()

Input: None

Output: Creates a global environment variable that contains the list of acceptable GIANT tissues
"""
def init_GIANT_queries():
    global GIANT_tissues
    GIANT_tissues = set()
    with open(PATH_TO_GIANT_TISSUES_LIST) as inFile:
        for line in inFile.readlines():
            GIANT_tissues.add(line.strip())

"""
get_GIANT_tissues_oi()

Input: 
iterable set of tissues of interest, as formatted by get_tissues_oi

Output: 
List of GIANT tissues that overlap with the tissues of interest
"""
def get_GIANT_tissues_oi(tissues):
    GIANT_tissues_oi = set()
    for tissue in tissues.keys():
        if tissue in GIANT_tissues:
            GIANT_tissues_oi.add(tissue)
    return GIANT_tissues_oi

"""
get_drug_tissues(drug_name)
Input:
drug_name - common drug name all lowercase

Output:
dictionary of tissues with the count of the number of 
"""
def get_drug_tissues(drug_name):
    ligand_targets = get_ligand_targets(drug_name)
    target_tissues = get_tissues_oi(ligand_targets)
    return(target_tissues)


# if __name__ == "__main__":
#     # get_drug_tissues("adapalene")
#     ligand_targets = get_ligand_targets("adapalene")
#     print(ligand_targets)
    # print(get_disease_id("acne vulgaris"))
    # qID = get_target_id(ligand_targets[0])
    # print(qID)

    # print(get_GO_terms('RARG'))
    # print(get_target_id('RARG'))
    # print(get_pathways('RARG'))
    # print(get_drug_tissues("adapalene"))
    # print(get_tissues_oi(["RARA"]))
    # temp = get_target_interactions(qID)
    # print("GENE TARGETS:", ligand_targets)
    # target_tissues = get_tissues_oi(ligand_targets)
    # print("GENE TISSUES:", target_tissues)
    # init_GIANT_queries()
    # GIANT_tiss = get_GIANT_tissues_oi(target_tissues)
    # print("GIANT TISSUES:", GIANT_tiss)


    """
    Things i want for iris:
    get_drug_tissues - takes in a common drug name, outputs a dictionary with tissue-count key value pairs
    get_tissues_oi - takes in a list of gene names, outputs a dictionary with tissue-count key value pairs
    """
