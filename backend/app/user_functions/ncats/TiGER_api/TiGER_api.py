import os
import pickle
import pandas as pd


TiGER_path = os.path.abspath(os.path.dirname(__file__))

TiGeR_reference_dir = os.path.join(TiGER_path, "TiGER_DB")

pkl_file = open(os.path.join(TiGeR_reference_dir, 'gene2Tissue_dict.pkl'), 'rb')
tissue_dict = pickle.load(pkl_file)
pkl_file.close()

def get_tissue_counts(geneList):
    results = dict()
    for gene in geneList:
        if gene in tissue_dict:
            for tissue in tissue_dict.get(gene):
                if tissue in results:
                    results[tissue] += 1
                else:
                    results[tissue] = 1
        else:
            if "all" in results:
                results["all"] += 1
            else:
                results["all"] = 1
    out = pd.DataFrame.from_dict(results, orient="index")
    out.columns = ["count"]
    out.sort_values(by = "count", inplace=True, ascending=False)
    return(out)



if __name__ == "__main__":
    test = get_tissue_counts(["ACE", "AGT", "REN"])
    test2 = get_tissue_counts(["AGT", "REN", "ACE","GNB3","WNK1","AGTR1","IGAN1","NF1","ATP2B1","FN1","NEDD4L","UTS2","UMOD","INS","EDN1","BTN2A1","LEP","INSR","CRP","ACE2","NOS3","KLK4","WNK4","PRSS8","HCCAT5","NPPA","CYP11B2","PTH","PEE1","CORIN"])
    print(test2)
