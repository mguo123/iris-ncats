import os
import pickle
import pandas as pd


TiGER_path = os.path.abspath(os.path.dirname(__file__))

ncats_path = os.path.dirname(os.path.dirname(TiGER_path))

TiGeR_reference_dir = os.path.join(ncats_path, "DB_data/TiGER_DB")

with open(os.path.join(TiGeR_reference_dir, 'gene2Tissue_dict.pkl'), 'rb') as f:
    tissue_dict = pickle.load(f)

with open(os.path.join(TiGeR_reference_dir, 'tissue_freq.pkl'), 'rb') as g:
    tissue_freq = pickle.load(g)
tissue_freq_values_avg = sum(tissue_freq.values()) / len(tissue_freq.values())


def get_tissue_counts(geneList):
    results = dict()
    total_genes = 0 # a sum of weighted gene numbers, genes are normalized based on frequency a tissue appears in the database
    for gene in geneList:
        if gene in tissue_dict:
            for tissue in tissue_dict.get(gene):
                if tissue in results:
                    results[tissue] += 1 / tissue_freq[tissue]
                    total_genes += 1 / tissue_freq[tissue]
                else:
                    results[tissue] = 1 / tissue_freq[tissue]
                    total_genes+=1 / tissue_freq[tissue]
        else:
            if "all" in results:
                results["all"] += 1 / tissue_freq_values_avg
                total_genes+=1 /tissue_freq_values_avg
            else:
                results["all"] = 1 /tissue_freq_values_avg
                total_genes+=1 /tissue_freq_values_avg

    # num_tissue_spec_genes = total_genes - all_genes

    # normalize tissue counts to the total of the tissue count then save the result as a string
    for tissue in results:
        results[tissue] =  "%.4f" % (results[tissue]/(float(total_genes) + 1e-10))
    
    out = pd.DataFrame.from_dict(results, orient="index")
    out.columns = ["Proportion of Genes"]
    out.index.name = 'Tissue'
    out.reset_index(inplace=True)

    out.sort_values(by ="Proportion of Genes", inplace=True, ascending=False)

    return(out)

if __name__ == "__main__":
    test = get_tissue_counts(["ACE", "AGT", "REN"])
    test2 = get_tissue_counts(["AGT", "REN", "ACE","GNB3","WNK1","AGTR1","IGAN1","NF1","ATP2B1","FN1","NEDD4L","UTS2","UMOD","INS","EDN1","BTN2A1","LEP","INSR","CRP","ACE2","NOS3","KLK4","WNK4","PRSS8","HCCAT5","NPPA","CYP11B2","PTH","PEE1","CORIN"])
    print(test2)
