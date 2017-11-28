import pickle

drugBank_targets = dict()

with open("/Users/alavertu/Desktop/Lab/NCATS_final/iris-ncats/backend/app/user_functions/ncats/shared_data/drugbank_targets.tsv") as inFile:
    i = 0
    inFile.readline()
    for line in inFile.readlines():
        info = line.strip().split("\t")
        if len(info) > 1:
            drugBank_targets[info[0].lower()] = info[1].split(",")
            i +=1
        # else:
        #     print(info)
        # if i > 10:
        #     break

with open('drugbank_target_dict.pkl', 'wb') as handle:
    pickle.dump(drugBank_targets, handle)
handle.close()
