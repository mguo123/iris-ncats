# globals.py
# loading data files
# print('loading graph files')

# global G_TO_RSIDS 
# G_TO_RSIDS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rsid_eQTL.pkl'),'rb'))
# global RS_TO_DATA 
# RS_TO_DATA = pickle.load(open(os.path.join(RSCS_DIR,'rsid_g_pval_rsqr.pkl'),'rb')) # rsid_g_pval_rsqr[rs]=[g,pv,rsq] 
# global G_TO_DISGENNET 
# G_TO_DISGENNET = pickle.load(open(os.path.join(RSCS_DIR,'disGeNet_gene_dis_score_dict.pkl'),'rb')) #
# global G_TO_OMIM
# G_TO_OMIM = pickle.load(open(os.path.join(RSCS_DIR,'OMIM_genes_to_phenotypes.pkl'),'rb')) # OMIM phenotypes
# global G_TO_PHW_SNPS 
# G_TO_PHW_SNPS = pickle.load(open(os.path.join(RSCS_DIR,'gene_to_rs_phWAS.pkl'),'rb')) # gene to SNPs from PheWAS
# global PHW_SNPS_TO_PHEN 
# PHW_SNPS_TO_PHEN = pickle.load(open(os.path.join(RSCS_DIR,'rs_to_phenOdds_phWAS.pkl'),'rb')) # SNPs to phenotype
# global ALL_ASSOC 
# ALL_ASSOC = pickle.load(open(os.path.join(RSCS_DIR,'all_assoc_to_nodes.pkl'),'rb'))# all associations, and their genes/SNPs
# global INTOME_SIZE 
# INTOME_SIZE = pickle.load(open(os.path.join(RSCS_DIR,'interactome_size.pkl'),'rb'))