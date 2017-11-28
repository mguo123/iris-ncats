# Written to get all associations with
# protein-protein network data
# written 3-23-17, JLW


#### USE settings.py to import all dependancies
# examples: from settings import *
import pickle,os
from collections import defaultdict
from scipy.stats import hypergeom


### RELEVANT PATHS ####
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # should get to the ncats directory
RSCS_DIR = os.path.join(PARENT_DIR, 'DB_data/PathFx_DB')
G_TO_RSIDS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rsid_eQTL.pkl'),'rb'))
RS_TO_DATA = pickle.load(open(os.path.join(RSCS_DIR,'rsid_g_pval_rsqr.pkl'),'rb')) # rsid_g_pval_rsqr[rs]=[g,pv,rsq] 
G_TO_DISGENNET = pickle.load(open(os.path.join(RSCS_DIR,'disGeNet_gene_dis_score_dict.pkl'),'rb')) #
G_TO_OMIM = pickle.load(open(os.path.join(RSCS_DIR,'OMIM_genes_to_phenotypes.pkl'),'rb')) # OMIM phenotypes
G_TO_PHW_SNPS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rs_phWAS.pkl'),'rb')) # gene to SNPs from PheWAS
PHW_SNPS_TO_PHEN = pickle.load(open(os.path.join(RSCS_DIR,'rs_to_phenOdds_phWAS.pkl'),'rb')) # SNPs to phenotype
ALL_ASSOC = pickle.load(open(os.path.join(RSCS_DIR,'all_assoc_to_nodes.pkl'),'rb'))# all associations, and their genes/SNPs
INTOME_SIZE = pickle.load(open(os.path.join(RSCS_DIR,'interactome_size.pkl'),'rb'))

def get_rsid_list(gene_list):
    genes_rsids = [[gene,rs] for gene in gene_list for rs in G_TO_RSIDS[gene]]
    return genes_rsids

def get_disease_associations(gene_list):
    dis_assoc = [[g,disease,dscore] for g in gene_list for (disease,dscore) in G_TO_DISGENNET[g].items() if g in G_TO_DISGENNET]
    return dis_assoc

def get_phenotype_associations(gene_list):
    pheno = []
    for g in gene_list:
        if g in G_TO_OMIM:
            if ';' in G_TO_OMIM[g]:
                for ph in G_TO_OMIM[g].split(';'):
                    pheno.append((g,ph))
            else:
                pheno.append((g,G_TO_OMIM[g]))
    return pheno

def get_pheWAS_SNPs(gene_list):
    genes_snps = list(set([(gene,snp) for gene in gene_list for snp in G_TO_PHW_SNPS[gene]]))
    return genes_snps

def get_snp_pheWAS(snp_list):
    snp_pheWAS = [[snp,PHW_SNPS_TO_PHEN[snp][0]] for snp in snp_list]
    return snp_pheWAS

def get_network_interactions(f):
    fdata = [l.strip().split('\t') for l in open(f,'r').readlines()]
    return fdata

def get_node_list(f):
    fdata = [l.strip().split('\t') for l in open(f,'r').readlines()]
    sourcen = [l[0] for l in fdata]
    sinkn = [l[1] for l in fdata]
    node_list = list(set(sourcen+sinkn))
    if '' in node_list:
        node_list.remove('')
    return node_list
    
def get_associations(netf, aname, rdir):
    '''
    Get all signficiant associations based on the hypergeometric probabilty of a path occuring

    inputs:
        <str> netf - path to network file
        <str> aname - unused
        <str> rdir - reference directory 

    outputs: none written to file
    
    '''

    print('gathering network data')
    # Gather network data
    node_list = get_node_list(rdir+netf)    
    net_int = get_network_interactions(rdir+netf)
    
    print('gathering associations')
    # pull assocations
    rsids = get_rsid_list(node_list) #eQTL
    disease_associations = sorted(get_disease_associations(node_list),key = lambda x:x[2],reverse=True)
    phenotypes = get_phenotype_associations(node_list)
    pheWAS_snps = get_pheWAS_SNPs(node_list)
    just_snps = [s for [g,s] in pheWAS_snps]
    snps_pheno = get_snp_pheWAS(just_snps)

    print('calculating hypergeometric probabilities')
    # Total up disease/phenotype associations, score
    # net_assoc = [d for [g,d,s] in disease_associations]+[d for [g,d] in phenotypes]+[p for [s,p] in snps_pheno] 
    assoc_count = defaultdict(int)
    for alist in disease_associations+phenotypes+snps_pheno:
        assoc_count[alist[1]]+=1
    N = INTOME_SIZE
    n = len(node_list)
    assoc_analy = []
    for (a,k) in assoc_count.items():
        K = len(ALL_ASSOC[a])
        prb = 1 - hypergeom.cdf(k,N,K,n)
        assoc_analy.append([a,k,K,prb])
    
    # Multiple hypothesis correct
    #Q = 0.05
    Q = 0.001
    sort_assoc = sorted(assoc_analy,key = lambda x:x[3])
    m = len(sort_assoc)
    mhc_assoc = []
    for (i,[a,k,K,prb]) in enumerate(sort_assoc):
        BH = (float(i+1)/m)*Q # calculate Benjamini-Hochberg based on ranked data
        mhc_assoc.append([i+1,a,k,K,prb,BH])
    sig_assoc = []
    rev_sort = sorted(mhc_assoc,key = lambda x: x[0],reverse=True)
    stop_rank = 0
    for [i,a,k,K,prb,BH] in rev_sort:
        if prb < BH:
            print('found stop-rank')
            stop_rank = i   
            break
    sig_assoc[:stop_rank]

    print('saving to output')


    # hypergeometric analysis
    network_name = os.path.basename(netf).split('.txt')[0]
    outf_name = os.path.join(rdir,network_name+'_assocations.txt')
    # outf = open(rdir+'/'+network_name+'_assoc_prb.txt','w')
    outf = open(outf_name,'w', encoding='utf-8')
    print('Total nodes: '+str(N))
    print('Neighborhood size: '+str(n))
    outf.write('\t'.join(['rank','phenotype','assoc in neigh','assoc in intom','probability','Benjamini-Hochberg','\n']))
    for [i,a,k,K,prb,BH] in mhc_assoc:
        new_line = '\t'.join([str(i),a ,str(k),str(K),str(prb),str(BH)+'\n'])
        outf.write(new_line)   
    outf.close()
