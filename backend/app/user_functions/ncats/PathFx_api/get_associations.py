# Written to get all associations with
# protein-protein network data
# written 3-23-17, JLW
# Modified margaret, 10/26/17


import pickle,os
# from optparse import OptionParser
from collections import defaultdict
from scipy.stats import hypergeom
from globals import *
# # import association file data
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # should yield the same as .. (with subfolders rcsc, scripts, and results)
RSCS_DIR = os.path.join(PARENT_DIR, 'shared_data')
G_TO_RSIDS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rsid_eQTL.pkl'),'rb'))
RS_TO_DATA = pickle.load(open(os.path.join(RSCS_DIR,'rsid_g_pval_rsqr.pkl'),'rb')) # rsid_g_pval_rsqr[rs]=[g,pv,rsq] 
G_TO_DISGENNET = pickle.load(open(os.path.join(RSCS_DIR,'disGeNet_gene_dis_score_dict.pkl'),'rb')) #
G_TO_OMIM = pickle.load(open(os.path.join(RSCS_DIR,'OMIM_genes_to_phenotypes.pkl'),'rb')) # OMIM phenotypes
G_TO_PHW_SNPS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rs_phWAS.pkl'),'rb')) # gene to SNPs from PheWAS
PHW_SNPS_TO_PHEN = pickle.load(open(os.path.join(RSCS_DIR,'rs_to_phenOdds_phWAS.pkl'),'rb')) # SNPs to phenotype
ALL_ASSOC = pickle.load(open(os.path.join(RSCS_DIR,'all_assoc_to_nodes.pkl'),'rb'))# all associations, and their genes/SNPs
INTOME_SIZE = pickle.load(open(os.path.join(RSCS_DIR,'interactome_size.pkl'),'rb'))
PHEN_TO_CUI = pickle.load(open(os.path.join(RSCS_DIR,'phenotype_to_cui.pkl',)'rb')) # phenotypes to cuis

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
    
def get_assoc(node_list):
    # pull assocationes
    #rsids = get_rsid_list(node_list) #eQTL
    disease_associations = sorted(get_disease_associations(node_list),key = lambda x:x[2],reverse=True) #[g,disease,score]
    phenotypes = get_phenotype_associations(node_list) # [(g,g_to_OMIM)]
    pheWAS_snps = get_pheWAS_SNPs(node_list)
    just_snps = [s for [g,s] in pheWAS_snps] 
    snps_pheno = get_snp_pheWAS(just_snps) # [snp,snpPheWAS]

    net_assoc = disease_associations+phenotypes+snps_pheno
    return net_assoc


def calc_hyp(net_assoc,all_assoc,N,Q,n):
    # count associations and store genes associated
    assoc_count = defaultdict(int)
    phen_genes = defaultdict(list)
    for alist in net_assoc:
        [phen,gene] = [alist[1],alist[0]]
        assoc_count[phen]+=1
        phen_genes[phen].append(gene)
    #n = len(node_list)
    assoc_analy = []
    for (a,k) in assoc_count.items():
        K = len(all_assoc[a])
        prb = 1 - hypergeom.cdf(k,N,K,n)    
        assoc_analy.append([a,k,K,prb])
    # Q = 0.001
    sort_assoc = sorted(assoc_analy,key = lambda x:x[3])
    m = len(sort_assoc)
    mhc_assoc = []
    for (i,[a,k,K,prb]) in enumerate(sort_assoc):
        BH = (float(i+1)/m)*Q # calculate Benjamini-Hochberg based on ranked data
        mhc_assoc.append([i+1,a,k,K,prb,BH])
    sig_assoc = []
    for [rank,phen,assnet,assint,prd,BH] in mhc_assoc:
        if prb<BH and assint >24:
        # if prb<BH:
            genes = phen_genes[phen]
            gene_str = ','.join(genes)
            if phen in PHEN_TO_CUI:
                cui = PHEN_TO_CUI[phen][0]
            else:
                cui = phen
            sig_assoc.append([rank,phen,cui,assnet,assint,prd,BH,gene_str])
        elif prb>BH:
            break
    return sig_assoc

def write_to_output(sig_assoc,outfname):
    outf = open(outfname,'w')
    outf.write('\t'.join(['rank','phenotype','cui','assoc in neigh','assoc in intom','probability','Benjamini-Hochberg','genes','\n']))
    if len(sig_assoc) >0:
        for line in sig_assoc:
            outf.write('\t'.join([str(x) for x in line])+'\n')
    outf.close()

    
def get_associations(netf, aname, rdir):
    # parser=OptionParser()

    # parser.add_option('-f','--file',dest='netf',help='Tab-separated network file of HUGO gene symbols')
    # parser.add_option('-a','--analysis_name',dest='aname',help='Name of analysis, will be appended to output files; experiment date is suggested')
    # parser.add_option('-d','--dir',dest='res_dir',help='Results directory. If none provided, a directory will be created matching the analysis name in the ../results/ dir')

    # (options,args) = parser.parse_args()

    # # check if results directory exists, otherwise assign based on analysis
    # if options.res_dir is None:
    #     rdir = '../results/'+options.aname+'/'
    # else:
    #     rdir = options.res_dir

    print('gathering network data')
    # Gather network data
    node_list = get_node_list(os.path.join(rdir,netf))
    net_int = get_network_interactions(os.path.join(rdir,netf))
    
    print('calculating hypergeometric probabilities')
    net_assoc = get_assoc(node_list)

    # Multiple hypothesis correct
    N = INTOME_SIZE
    Q = 0.001
    n = len(node_list)
    sig_assoc = calc_hyp(net_assoc,ALL_ASSOC,N,Q,n) 
    network_name = os.path.basename(netf).split('.txt')[0]
    outfname = os.path.join(rdir,network_name+'_assocations.txt')
    print ('saving to output', outfname)

    write_to_output(sig_assoc,outfname)

    return sig_assoc
