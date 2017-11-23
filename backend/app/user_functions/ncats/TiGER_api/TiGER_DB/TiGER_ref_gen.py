import os
import pickle

TiGER_path = os.path.abspath(os.path.dirname(__file__))

i = 0



################# GET ID TO HS MAPPINGS #################
gene2hs_db = dict()

with open(os.path.join(TiGER_path, "symbol2hs-Table.txt")) as gene_reference:
    gene_reference.readline()
    for line in gene_reference.readlines():
        info = line.strip().split("\t")
        for id in info[1:]:
            gene2hs_db[id] = info[0]
        # i +=1

# print(gene2hs_db)

# with open('gene2hs_dict.pkl', 'wb') as handle:
#     pickle.dump(gene2hs_db, handle)
# handle.close()


################# GET TISSUE EXP DATA #################
TiGER_tissue_DB = dict()
with open(os.path.join(TiGER_path, "hs2tissue-Table.txt")) as gene_reference:
    gene_reference.readline()
    for line in gene_reference.readlines():
        info = line.strip().split("\t")

        gene = gene2hs_db.get(info[0])
        if gene is not None:
            if gene in TiGER_tissue_DB:
                [TiGER_tissue_DB[gene].add(x) for x in info[1:]]
            else:
                TiGER_tissue_DB[gene] = set(info[1:])
        # i +=1
        # if i > 10:
        #     break

with open('gene2Tissue_dict.pkl', 'wb') as handle:
    pickle.dump(TiGER_tissue_DB, handle)
handle.close()
