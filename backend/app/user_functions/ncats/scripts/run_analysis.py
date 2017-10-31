# written to contain some of the core
# methods for analyzing neighborhoords and post-processing
# written 7-13-17, JLW
# Modified margaret, 10/26/17


import csv, pickle, os, sys, find_neighborhood_beta
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import pandas as pd
from textwrap import wrap
from collections import defaultdict
from find_neighborhood_beta import find_neighborhood as fgn
# import get_associations
from get_associations_deprecated import get_associations

def write_neighborhood_to_file(pth_dic,outf):
    for (pth,pscore) in pth_dic.items():
        if '@' in pth:
            a=pth.split('@')[-2]
            b=pth.split('@')[-1]
            outf.write('\t'.join([a,b,str(pscore),'\n']))
        else:
            a=pth
            outf.write('\t'.join([a,'',str(pscore),'\n']))
    outf.close()

# need to set rand_dir 
def calculate_specificity(pth_dic,outf,outf2,scr_thr):
    gene_scores = dict([(pth.split('@')[-1],pscore) for (pth,pscore) in pth_dic.items()])
    neigh_genes = [k for k in gene_scores.keys()]
    enriched_genes = []
    outf.write('\t'.join(['gene','score list','average','enrichment','\n']))
    for tgene in neigh_genes:
        # update tgene for weird text-isms
        sn  = tgene # the save-name for the gene
        if '/' in sn or '+' in sn or ':' in sn:
            sn = sn.replace('+','-')
            sn = sn.replace(':','-')
            sn = sn.replace('/','-')
        grf = rand_dir+'_'.join([sn,str(scr_thr),'randPathScores','.pkl'])
        if not os.path.exists(grf): # genes that have no specificity scores (i.e. no connections above the threshold)
            continue
        gpathlist = pickle.load(open(grf,'rb'))
        # average
        avg = np.mean(gpathlist)
        # enrichment
        enr = gene_scores[tgene]-avg

        outf.write('\t'.join([tgene,str(avg),str(enr),'\n']))
        if enr >0:
            enriched_genes.append(tgene)
    outf.close()

    spec_pth_dic = {}
    for (pth,pscore) in pth_dic.items():
        if '@' in pth:
            b=pth.split('@')[-1]
            if b in enriched_genes:
                spec_pth_dic[pth]=pscore
    write_neighborhood_to_file(spec_pth_dic,outf2)

    return enriched_genes

def merge_networks(allf,rdir,drug):
    # merge networks
    print('merging protein networks')
    nfs = [f for f in allf if 'specific_neighborhood.txt' in f]
    outfname = drug+'_merged_neighborhood.txt'
    outf = open(os.path.join(rdir,outfname),'w')
    unique_edges = set()
    for n in nfs:
        flines = [l.strip() for l in open(rdir+n,'rU').readlines()]
        for l in flines:
            a =l.strip().split('\t')[0]
            b = l.strip().split('\t')[1]
            if (a,b) not in unique_edges or (b,a) not in unique_edges:
                n = outf.write(l)
                n = outf.write('\n')
                unique_edges.add((a,b))
    outf.close()
    return outfname

def check_if_drug_in_network(dts, netxobj):
    print('Checking targets')
    all_drug_targets = np.array([])
    not_found = np.array([])    
    for (drug,target_list) in dts:
        # print(drug)
        # print(target_list)
        all_drug_targets = np.append(all_drug_targets, target_list)
        # print(netxobj.__dict__.keys())
        # print(find_neighborhood_beta)
        for i, drug_target in enumerate(target_list):
            # if drug_target not in netxobj.nodes_iter():
            if drug_target not in find_neighborhood_beta.G:
                print(drug_target, 'not found in network of drug', drug)
                not_found = np.append(not_found, drug_target)
    all_drug_targets = np.unique(all_drug_targets)
    not_found = np.unique(not_found)
    print(not_found.shape[0], 'not found of ', all_drug_targets.shape[0], 'total targets')
    
def run_all_drugs(dts,rdir,netxobj,scr_thr):
    find_neighborhood_beta.G = netxobj
    all_merge_files = []
    # print('dts', dts)
    for (drug,target_list) in dts:
        print(drug, target_list)
        # print(netxobj.__dict__.keys())
        # print(find_neighborhood_beta)
        for drug_target in target_list:
            # if drug_target not in netxobj.nodes_iter():
            if drug_target not in find_neighborhood_beta.G:
                print(drug_target, 'not found in network')
                continue        
            print('Gene: '+drug_target)
            res_dir = os.path.join(rdir,drug+'_networks/')
            if not os.path.exists(res_dir):
                os.mkdir(res_dir)
            
            # create network
            pth_dic = fgn(drug_target,scr_thr)
            
            # save network
            with open(os.path.join(res_dir, drug_target+str(scr_thr)+'_neighborhood.txt'),'w') as outf:
                write_neighborhood_to_file(pth_dic,outf)
            
            if len(pth_dic) > 1: # for target networks that extend beyond the target
                # calculate specificity
                outfname = os.path.join(res_dir, drug_target+'_thr_'+str(scr_thr)+'_specificity.txt')
                with open(outfname,'w') as outf:

                    outf2name =os.path.join(res_dir, drug_target+str(scr_thr)+'_specific_neighborhood.txt')
                    with open(outf2name,'w') as outf2:
                        enriched_genes = calculate_specificity(pth_dic,outf,outf2,scr_thr)
            else: # for target networks that only contain the target (i.e. the target was not interacting with any proteins above the threshold)
                outf2name = os.path.join(res_dir, drug_target+str(scr_thr)+'_specific_neighborhood.txt')
                with open(outf2name,'w') as outf2:
                    write_neighborhood_to_file(pth_dic,outf2)
        
        # merge files, do phenotype enrichment
        # print('here')
        allf = [f for f in os.listdir(res_dir)]
        # print(allf)
        spnn = merge_networks(allf,res_dir,drug)
        # print('merged networks')
        aname = 'merged'
        sig_assoc = get_associations(spnn, aname, res_dir)

        all_merge_files.append(res_dir+spnn)

    return all_merge_files

# method for ranking rows
def msum(l):
    newl = l[np.nonzero(l)]
    return np.mean(newl)

# make a heatmap
def make_a_heatmap(drugs,all_merge_files,rdir,analysis_name,per_cutoff):
    data_dic = defaultdict(list)
    ndrugs = len(drugs)
    dind = dict([(d,i) for (i,d) in enumerate(drugs)])
    for mf in all_merge_files:
            (pth,fname) = os.path.split(mf)
            mpf = os.path.join(pth,fname.split('.txt')[0]+'_assoc_prb.txt')
            drug = fname.split('_')[0]
            dR = csv.DictReader(open(mpf,'rU'),delimiter='\t')
            for l in dR:
                    if float(l['probability']) < float(l['Benjamini-Hochberg']):
                            if float(l['assoc in intom']) > 24:
                                    phen = l['phenotype']
                                    asn = float(l['assoc in neigh'])
                                    asi = l['assoc in intom']
                                    pkey = ' '.join([phen,'('+asi+')']) # key for the heatmap
                                    pern = asn/float(asi) # percent of network
                                    i = dind[drug]

                                    if pern > per_cutoff:
                                            if pkey in data_dic:
                                                    data_dic[pkey][(0,i)] = pern
                                            else:
                                                    zarray = np.zeros((1,ndrugs))
                                                    zarray[(0,i)] = pern
                                                    data_dic[pkey] = zarray
                    else:
                            break

    pickle.dump(data_dic,open(rdir+'merged_phen_dic.pkl','wb'))

    # format into data frame
    columns = drugs
    rank_data = sorted([(k,msum(v)) for (k,v) in data_dic.items()],key = lambda x:x[1],reverse=True)
    rank_phen = [k for (k,s) in rank_data]
    print('how many phenotypes above the threshold?')
    print(len(rank_phen))
    new_array = data_dic[rank_phen[0]]
    for k in rank_phen[1:]:
        new_array = np.concatenate((new_array,data_dic[k]),axis=0)

    df1 = pd.DataFrame(new_array,index=rank_phen,columns=columns)
    vmax = np.amax(new_array)
    vmin = np.amin(new_array)

    # plot it
    print('plotting')
    fix,ax = plt.subplots()
    heatmap = ax.pcolor(df1,cmap=plt.cm.plasma,alpha=0.8,vmin=vmin,vmax=vmax)
    fig = plt.gcf()
    fig.set_size_inches(8,11)

    # format
    ax.set_frame_on(False)
    ax.set_yticks(np.arange(df1.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(df1.shape[1]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    row_labels = ['\n'.join(wrap(l,70)) for l in df1.index]
    ax.set_xticklabels(columns, minor=False,fontsize=9)
    ax.set_yticklabels(row_labels, minor=False,fontsize=7.5)
    plt.xticks(rotation=90)
    ax.grid(False)

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    plt.gcf().subplots_adjust(left=0.6)
    plt.gcf().subplots_adjust(top=0.65)
    cbar = plt.colorbar(heatmap)
    cbar.ax.tick_params(labelsize=9)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('percent engaged', rotation=270,fontsize=9)

    plt.savefig(rdir+'over30_ranked_'+analysis_name+'_phenotypes.png',format='png')

