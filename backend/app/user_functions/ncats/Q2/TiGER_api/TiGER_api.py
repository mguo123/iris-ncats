import os
import pickle
import pandas as pd


TiGER_path = os.path.abspath(os.path.dirname(__file__))

ncats_path = os.path.dirname(os.path.dirname(TiGER_path))

TiGeR_reference_dir = os.path.join(ncats_path, "DB_data/TiGER_DB")

pkl_file = open(os.path.join(TiGeR_reference_dir, 'gene2Tissue_dict.pkl'), 'rb')
# print('TIGER!!!!!')
# print(ncats_path)
# print(pkl_file)

tissue_dict = pickle.load(pkl_file)
pkl_file.close()

def get_tissue_counts(geneList):
    results = dict()
    total_genes = 0
    all_genes = 0
    for gene in geneList:
        if gene in tissue_dict:
            for tissue in tissue_dict.get(gene):
                if tissue in results:
                    results[tissue] += 1
                    total_genes+=1
                else:
                    results[tissue] = 1
                    total_genes+=1
        else:
            if "all" in results:
                results["all"] += 1
                total_genes+=1
                all_genes+=1
            else:
                results["all"] = 1
                total_genes+=1
                all_genes+=1

    # num_tissue_spec_genes = total_genes - all_genes

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
