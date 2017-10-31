from iris import state_types as t
from iris import IrisCommand


from iris import state_machine as sm
from iris import util as util
from iris import iris_objects

class DrugDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "find the mechanism of action of {drug} for treating {disease}"
    
    # give an example for iris to recognize the command
    examples = ["mechanism of action", "treat disease", "how drug works", "how drug works on disease", "how does drug affect disease", "how does {drug} affect {disease}" "How does {drug} treat {disease}", "What is the mechanism of action for {drug} treating {disease}", "is {disease} treatable by {drug}"]
    
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"drug":t.String("What is the drug you want to analyze?"), "disease":t.String("What is the disease you want to analyze?")}
    
    # core logic of the command
    def command(self, drug, disease):
        # import numpy
        # return numpy.random.randint(100)
        import sys
        print (sys.path)
        sys.path.append("/Users/margaret/Documents/iris-agent") #### NEED TO FIX!!!! 
        print (sys.path)
        from node_modules.ncats.scripts import run_test
        return run_test.run_drug_single(drug, disease)
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        # return ["Book a flight to DC to get the answer. Here's your lucky number", result]
        return ["Answer was found: " + result[0], "Answer: ", result[1]]

_DrugDisease = DrugDisease()


class DrugDiseaseMulti(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "find the mechanism of action of this list of drug-disease pairs"
    
    # give an example for iris to recognize the command
    examples = ["multiple drugs and diseases", "how list of drugs works", "how multiple drug works on disease", "how do these drugs affect these diseases", "how does this list of {drug_list} affect {disease_list}" ]

    # type annotations for each command argument, to help Iris cosllect missing values from a user
    argument_types = {"drug_list":t.List("What drugs do you want to analyze?"), "disease_list":t.List("What diseases do you want to analyze? List must be the same length as the list of drugs.")}
    
    # core logic of the command
    def command(self, drug_list, disease_list):
        # import numpy
        # return numpy.random.randint(100)
        assert(len(drug_list) == len(disease_list))
        import sys
        sys.path.append("/Users/margaret/Documents/iris-agent") #### NEED TO FIX!!!! 
        from node_modules.ncats.scripts import run_test
        results = []
        for drug, disease in zip(drug_list, disease_list):
            result = run_test.run_drug_single(drug, disease)
            result_txt = ' '.join(('Answer Found:', result[0], 'Answer:', result[1]))
            results.append(' '.join((drug, disease, result_txt)))
        return '\t'.join(results)
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, results):

        # return ["Book a flight to DC to get the answer. Here's your lucky number", result]

        return results.split('\t')

_DrugDiseaseMulti = DrugDiseaseMulti()



class DrugDiseaseCSV(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "drug-disease pairs from csv {paired_csv} and saved to {saved_dir}"
    
    # give an example for iris to recognize the command
    examples = ["multiple drugs and diseases from csv", "csv of drug disease matches", "loaded drug disease match", "find the mechanism of action of drug-disease pairs from csv {paired_csv} and saved to {saved_dir}"]

    # type annotations for each command argument, to help Iris cosllect missing values from a user
    argument_types = {"paired_csv":t.File("What drugs do you want to analyze?"), "saved_dir":t.String("Where do you want to save the results")}

    # core logic of the command
    def command(self, paired_csv, saved_dir):
        import pandas as pd
        new_df = iris_objects.IrisDataframe(None, empty=True)
        new_df.df = pd.read_csv(paired_csv.path, sep='\t')
        # print('read in file')
        self.iris.add_to_env('drug_disease_list', new_df)
        # print('added to environment')

        import sys
        sys.path.append("/Users/margaret/Documents/iris-agent") #### NEED TO FIX!!!! 

        from node_modules.ncats.scripts import run_test
        print('path', paired_csv.path)
        print('storage_dir', saved_dir)
        run_test.run_drug_multi(paired_csv.path, storage_dir = saved_dir, iterate_until_found=False)
        return ' '.join(("written results to:", saved_dir))
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, results):

        return results

_DrugDiseaseCSV = DrugDiseaseCSV()


