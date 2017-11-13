

# This function attempts to answer the question:
# Do any genetic diseases protect against this condition
def Q1_query(condition, genetic_disease=None):
    if condition.lower() == "malaria":
        results = ['Sickle cell disease protects against malaria']
        results.append("No evidence for association exists at the moment")
    elif condition.lower() == "ebola":
        results = ['Niemann-Pick Type C protects against ebola']
        results.append("No evidence for association exists at the moment")
    elif condition.lower() == "cholera":
        results = ['Cystic Fibrosis protects against ebola']
        results.append("No evidence for association exists at the moment")
    else:
        results = ["It is unknown whether any genetic disease protect against %s" % condition]
    return results