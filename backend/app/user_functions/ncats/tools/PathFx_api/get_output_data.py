# written to get all genes
# associated with a given phenotype
# for the gold standard drugs
# writen 8-17-17 JLW
# Modified margaret, 10/26/17


import csv, os, pickle
from collections import defaultdict

# # resources for pulling associations
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # should get to the ncats directory
RSCS_DIR = os.path.join(PARENT_DIR, 'DB_data/PathFx_DB')
g_to_disGenNet = pickle.load(open(RSCS_DIR+'/disGeNet_gene_dis_score_dict.pkl','rb'))
g_to_OMIM = pickle.load(open(RSCS_DIR+'/OMIM_genes_to_phenotypes.pkl','rb')) # OMIM phenotypes
g_to_phW_snps = pickle.load(open(RSCS_DIR+'/gene_to_rs_phWAS.pkl','rb')) # gene to SNPs from PheWAS
phW_snps_to_phen = pickle.load(open(RSCS_DIR+'/rs_to_phenOdds_phWAS.pkl','rb')) # SNPs to phenotype

def get_nodes(f):
    d = [l.strip().split('\t') for l in open(f,'rU').readlines()]
    a = [x[0] for x in d]
    b = [x[1] for x in d]
    return list(set(a+b))

def get_disease_associations(gene_list):
    dis_assoc = [[g,disease,dscore] for g in gene_list for (disease,dscore) in g_to_disGenNet[g].items() if g in g_to_disGenNet]
    return dis_assoc

def get_phenotype_associations(gene_list):
    pheno = []
    for g in gene_list:
        if g in g_to_OMIM:
            if ';' in g_to_OMIM[g]:
                for ph in g_to_OMIM[g].split(';'):
                    pheno.append((g,ph))
            else:
                pheno.append((g,g_to_OMIM[g]))
    return pheno

def get_pheWAS_SNPs(gene_list):
    genes_snps = list(set([(gene,snp) for gene in gene_list for snp in g_to_phW_snps[gene]]))
    return genes_snps

def get_snp_pheWAS(snp_list):
    snp_pheWAS = [[snp,phW_snps_to_phen[snp][0]] for snp in snp_list]
    return snp_pheWAS

def get_output_data(storage_dir, drug):
    """
    Get the output dictionary with genes relevant for each phenotype that is drug-relevant

    Input:
        <str> storage_dir: path
        <str> drug

    Output:
        <dictionary> assoc_to_genes - key: phenotype: value: list of genes that are associated with phenotype

    """
    drug_dir = os.path.join(storage_dir, drug + "_networks")

    drug_file = os.path.join(drug_dir, drug + '_merged_neighborhood.txt')
    print('drug file', drug_file)
    allf = [(d,sd,sf) for (d,sd,sf) in os.walk(storage_dir)] # dir, sub-dirs, dir_files

    merfs = [(sf.split('_')[0],d,sf) for (d,sd,sflist) in allf for sf in sflist if 'merged_neighborhood.txt' in sf] # just the merged neighborhood files
    if os.path.isfile(drug_file): 
        n = get_nodes(drug_file)
        disease_associations = sorted(get_disease_associations(n),key = lambda x:x[2],reverse=True)    
        phenotypes = get_phenotype_associations(n)
        pheWAS_snps = get_pheWAS_SNPs(n)
        just_snps = [s for [g,s] in pheWAS_snps]
        snps_pheno = get_snp_pheWAS(just_snps)

        assoc_to_genes = defaultdict(list)
        for phendata in disease_associations + phenotypes + snps_pheno:
            (gene,phenotype) = [phendata[0],phendata[1]] # just unpack the first two elements
            assoc_to_genes[phenotype].append(gene)

        assoc_to_genes.pop('',None)
        outf = os.path.join(drug_dir,drug+'_phens_to_genes.pkl')
        print('creating file: ', outf)
        pickle.dump(assoc_to_genes,open(outf,'wb'))
        return assoc_to_genes


