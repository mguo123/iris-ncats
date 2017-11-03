# written to create a better summary table
# from phenotypic enrichment, written
# 8-31-17 JLW
# Modified margaret, 10/26/17


import csv, os, pickle
# from optparse import OptionParser
from collections import defaultdict


def get_sum(mf,gtphen): # pass a merged phenotype result file
    sig_phens_data = []
    gtpd = pickle.load(open(gtphen,'rb'))
    dR = csv.DictReader(open(mf,'rU'),delimiter='\t')
    for l in dR:
        # gather relevant data
        [ph,rank,BH,prb] = [l['phenotype'],l['rank'],l['Benjamini-Hochberg'],l['probability']]
        [asii,asin] = [l['assoc in intom'],l['assoc in neigh']]    
        pern = float(asin)/float(asii)

        if float(prb) < float(BH) and float(asii)>10:
            sig_genes = gtpd[ph]
            sig_phens_data.append((ph,rank,BH,prb, asii,asin,str(pern),sig_genes))

        elif float(prb) > float(BH):
            break # exit the loop as soon as stop criteria are met
    return sig_phens_data


def get_results(storage_dir, drug, old_format = False, save_file = True):
    # parser=OptionParser()

    # parser.add_option('-d','--rdir',dest='rdir',help='The results directory where all drugs are sub-directories')
    # (options,args) = parser.parse_args()

    # get relevant phenotype summary and phens_to_genes dic
    # rdir = options.rdir
    # print(rdir)
    # allf = [(d,sd,sf) for (d,sd,sf) in os.walk(rdir)] # dir, sub-dirs, dir_files
    # # print(allf)
    # # merfs = [(sf.split('_')[0],d,sf) for (d,sd,sflist) in allf for sf in sflist if 'merged_neighborhood_assoc_prb.txt' in sf] # just the merged neighborhood files
    # merfs = [(sf.split('_')[0],d,sf) for (d,sd,sflist) in allf for sf in sflist if 'merged_neighborhood_assocations.txt' in sf] # just the merged neighborhood files
    
    drug_dir = os.path.join(storage_dir, drug + "_networks")
    if old_format:
        drug_file = os.path.join(drug_dir, drug + '_merged_neighborhood_assoc_prb.txt')

    else:
        drug_file = os.path.join(drug_dir, drug + '_merged_neighborhood_assocations.txt')
    print(drug)
    if os.path.isfile(drug_file): 
    # for (drug,drug_dir,drug_file) in merfs:
        # print(drug)
        mapf = drug_file
        gtphen = os.path.join(drug_dir,drug+'_phens_to_genes.pkl')
        sum_asscs = get_sum(mapf,gtphen)

        outfname = os.path.join(drug_dir,drug+'_merged_assc_full_summary.txt')
        
        if save_file:
            outf = open(outfname,'w')

            print('writing output to', outfname)
            outf.write('\t'.join(['phenotype','rank','BHcorrPval', 'Pval', 'assoc_in_intom','assoc_in_neigh','perc_overlap','neigh_genes_in_phen\n']))
            # for line in sum_asscs:
            #     print(line)
            for line in sum_asscs:
                if line is not None:
                    # print(len(line), line[0])
                    [ph,rank,BH,prb, asii,asin,pern,sig_genes] = line
                    sig_gn_str = ','.join(sig_genes)
                    outf.write('\t'.join([ph,rank,BH,prb, asii,asin,pern,sig_gn_str,'\n']))        
            outf.close()
        # print(sum_asscs)
        return sum_asscs

    else:
        return None


# if __name__ == "__main__":
#         main()
