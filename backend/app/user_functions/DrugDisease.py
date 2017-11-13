from iris import state_types as t
from iris import IrisCommand


from iris import state_machine as sm
from iris import util as util
from iris import iris_objects

from app.user_functions.ncats.scripts import run_test, run_main

import numpy as np
# run_test.run_drug_single('a', 'b')

class DrugDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "How does {drug} treat {disease}?"
    
    # give an example for iris to recognize the command
    examples = ["why does {drug} work against {disease}", "why does {drug} treat disease", "what are the targets of {drug} relevant to {disease}", "mechanism of action", "treat disease", "how drug works", "how drug works on disease", "how does drug affect disease", "how does {drug} affect {disease}" "How does {drug} treat {disease}", "What is the mechanism of action for {drug} treating {disease}", "is {disease} treatable by {drug}"]
    
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"drug":t.String("Just for confirmation: What is the drug you want to analyze?"), "disease":t.String("What is the disease you want to analyze?")}
    
    # core logic of the command
    def command(self, drug, disease):
        # import numpyfrom ncats.scripts import run_te
        # return numpy.random.randint(100)

        answer= run_main.run_drug_single(drug, disease)
        # print (answer)
        # print('!!!!!!!!!')
        return answer
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        # return ["Book a flight to DC to get the answer. Here's your lucky number", result]
        answer_found_bool, answer_found_exp_str, ph_genes_str, drug = result
        if answer_found_bool:
            # FROM backend, ph_genes_str = '\t'.join([prb, BH, disease, sig_genes]) # was done this way to save results in csv (for multi run)
            [prb, BH, disease, sig_genes] = ph_genes_str.split('\t')
            if float(prb) < float(BH):
                prob_string = ' '.join([drug, 'was found to treat', disease, 'with probability:', prb, 'for Benjamini–Hochberg significance level', BH])
                sig_genes_string = ' '.join(['Signficant genes were found as:', sig_genes])
            else:
                prob_string = ' '.join([drug, 'was found to treat', disease, ', but probability', prb, 'was not significant for Benjamini–Hochberg singificance level', BH])
                sig_genes_string = ' '.join(['Associated genes were found as:', sig_genes])

            return ["Answer was found!", prob_string, sig_genes_string, 'find the results in results folder']
        else:
            answer_explanation_string = [' '.join(['Answer was not found because:', answer_found_exp_str])]
            multi_answer_line = ['We queried the gene neighborhood of drug targets and found the following phenotypes to be significant', 'Here we list significant phenotypes in order of probability. Column headings are probability, significance level cutoff, phenotype, and a list of genes associated']
            result = answer_explanation_string + multi_answer_line

            ph_genes_arr = ph_genes_str.split('\t') # prb, BH, ph, sig_genes 
            ph_genes_array_grouped = [ph_genes_arr[x:x+4] for x in range(0, len(ph_genes_arr),4)]
            ph_genes_array_grouped_iris = iris_objects.IrisDataframe(column_names=["probability", "Benjamin Hochberg significance cutoff", "Phenotype", "list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_grouped)
            self.iris.add_to_env('drug_disease_results', ph_genes_array_grouped_iris)
            result.append(ph_genes_array_grouped_iris)  

            return result

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

        from ncats.scripts import run_test
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

