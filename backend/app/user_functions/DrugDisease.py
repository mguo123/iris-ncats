from iris import state_types as t
from iris import IrisCommand


from iris import state_machine as sm
from iris import util as util
from iris import iris_objects

from app.user_functions.ncats.Q2.PathFx_api import run_main
from app.user_functions.ncats.Q2.Q2_main import Q2_query
import numpy as np
import os
# run_test.run_drug_single('a', 'b')

overall_path = os.path.abspath(os.path.dirname(__file__))
results_dir = os.path.join(overall_path, "ncats/results/Q2/")

class Options(object):
    def __init__(self, verbose=True, batchFile=None, gen_image=False, gen_pubmed=False, outPath=results_dir):
        self.verbose= verbose
        self.gen_image = gen_image
        self.gen_pubmed = gen_pubmed
        self.outPath = outPath
        self.batchFile = batchFile


class DrugDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "How does {drug} treat {disease}?"
    
    # give an example for iris to recognize the command
    examples = ["why does {drug} work against {disease}", "why does {drug} treat {disease}", "what are the targets of {drug} relevant to {disease}", "mechanism of action", "treat {disease}", "how drug works", "how drug works on disease", "how does {drug} affect {disease}", "how does {drug} affect {disease}" "How does {drug} treat {disease}", "What is the mechanism of action for {drug} treating {disease}", "is {disease} treatable by {drug}"]
    
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"drug":t.String("Okay, a couple more questions to set up this task. For confirmation: What is the drug you want to analyze?"), 
                        "disease":t.String("What is the disease you want to analyze?"),
                        "bool_image":t.YesNo("Would you like to visualize the results as a diagram?",
                                    yes=True, no=False),
                        "bool_pubmed":t.YesNo("Would you like to get the list of pubmed IDs for reference?",
                                    yes=True, no=False),
                        "bool_other_disease":t.YesNo("Would you like to know other diseases that can be treated by this drug?",
                                    yes=True, no=False)
                            
                        }
    
    # core logic of the command
    def command(self, drug, disease, bool_image, bool_other_disease, bool_pubmed):
        print('BEFORE QUERY!!!')
        # create options structure 

        options = Options(gen_image=bool_image, gen_pubmed=bool_pubmed)


        answer = Q2_query(drug, disease, options)

        print('answered Q2 query!!!')
        
        if isinstance(answer, str):
            return answer

        if bool_other_disease:
            answer["other_disease"] = run_main.find_drug_indications(drug)

        answer['drug'] = drug
        answer['disease'] = disease
    
        return answer
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        # Components of result are in dictionary form:
        # result = {"GOENRICH":result, "drug_genes":drug_genes, "disease_genes":dis_genes, "dis_tissue_data":tissue_df, "dis_tissue_data_short":tissue_df_short, "image_file": image_path, "other_disease": jenn's result, "pubmed": PMIDs
        
        if isinstance(result, str):
            return result

        query_name = result['drug'][:min(len(result['drug']), 3)] + "_" + result['disease'][:min(len(result['disease']), 3)]
        query_name = '_'.join(query_name.split(' '))
        query_name = "_" + query_name
        print(query_name, "== query name")

        # Print out genes associated with drug
        result_array = []
        result_array.append('Top genes found to be targetted by %s are below. Full dataset saved as drug_genes_{drug_disease}' % result['drug'])
        drug_gene_term_object = iris_objects.IrisDataframe(data=result['drug_genes'])
        self.iris.add_to_env('drug_genes' + query_name, drug_gene_term_object) 
        drug_gene_term_object_short = iris_objects.IrisDataframe(data=result['drug_genes_short'])
        result_array.append(drug_gene_term_object_short)
        # result_array.append("Full dataset saved as drug_associated_genes")


        # Print out genes associated with disease
        result_array.append('Top genes found to be associated with %s are below. Full dataset saved as disease_genes_{drug_disease}' % result['disease'])        
        disease_gene_term_object = iris_objects.IrisDataframe(data=result['disease_genes'])
        self.iris.add_to_env('disease_genes' + query_name, disease_gene_term_object)
        disease_gene_term_object_short = iris_objects.IrisDataframe(data=result['disease_genes_short'])
        result_array.append(disease_gene_term_object_short)


        # Print out signficant go terms
        result_array.append('Top significant GO terms associated with the drug-disease interaction are shown. Full dataset saved as go_terms_{drug_disease}')
        go_term_object = iris_objects.IrisDataframe(data=result['GOENRICH'])
        self.iris.add_to_env('go_terms' + query_name, go_term_object)
        go_term_object_short = iris_objects.IrisDataframe(data=result['GOENRICH_short'])
        result_array.append(go_term_object_short)
        # result_array.append("Full dataset saved as drug_disease_go_terms")


        # get tissue = disease
        result_array.append('Top disease-relevant tissues in which this drug-disease interaction takes place are shown. Full dataset saved as tissues_{drug_disease} ')
        tissue_object_dis = iris_objects.IrisDataframe(data=result['tissue_df_dis'])
        tissue_object_dis_short = iris_objects.IrisDataframe(data=result['tissue_df_dis_short'])
        self.iris.add_to_env('tissues_disease' + query_name, tissue_object_dis)
        result_array.append(tissue_object_dis_short)
 
         # get tissue = drug
        result_array.append('Top drug-relevant tissues in which this drug-disease interaction takes place are shown. Full dataset saved as tissues_{drug_disease} ')
        tissue_object_drug = iris_objects.IrisDataframe(data=result['tissue_df_drug'])
        tissue_object_drug_short = iris_objects.IrisDataframe(data=result['tissue_df_drug_short'])
        self.iris.add_to_env('tissues_drug' + query_name, tissue_object_drug)
        result_array.append(tissue_object_drug_short)


        # display image
        if "image_file" in result:
            os.system("open " + result["image_file"])


        if "pubmed" in result:
            if isinstance(result["pubmed"], str):
                result_array.append(result["pubmed"])
            else:
                result_array.append("Following are PMIDs that support the interaction: Full dataset saved as pmid_{drug_disease}.")        
                pmid_df_short = iris_objects.IrisDataframe(data=result["pubmed_short"])
                pmid_df = iris_objects.IrisDataframe(data=result["pubmed"])
                self.iris.add_to_env('pmid' + query_name, pmid_df)
               
                result_array.append(pmid_df_short)
                # result_array.append("Full dataset saved as pmid_ids")

        # get other possible disease 
        if "other_disease" in result:
            ph_genes_str, drug = result["other_disease"]
            multi_answer_line = ['Top hits of diseases potentially impacted by %s. Full dataset saved as drug_indications_{drug_disease}.' % result['drug'], 'We queried the gene neighborhood of drug targets and found the following phenotypes to be significant. Here we list significant phenotypes in order of probability. Column headings are phenotype, probability, significance level cutoff, and a list of genes that support the relationship']
            result_array = result_array + multi_answer_line
            ph_genes_arr = ph_genes_str.split('\t') # prb, BH, ph, sig_genes 
            ph_genes_array_all = [ph_genes_arr[x:x+4] for x in range(0, len(ph_genes_arr),4)]
            ph_genes_array_all_iris = iris_objects.IrisDataframe(column_names=[ "Phenotype", "probability", "Benjamin Hochberg significance cutoff","list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_all)
            self.iris.add_to_env('drug_indications' + query_name, ph_genes_array_all_iris)
            ph_genes_array_short = [ph_genes_arr[x:x+4] for x in range(0, min(5*4,len(ph_genes_arr)),4)]
            ph_genes_array_short_iris = iris_objects.IrisDataframe(column_names=["Phenotype", "Probability", "Benjamin Hochberg significance cutoff", "list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_short)
            result_array.append(ph_genes_array_short_iris)  
            # result_array.append("Full dataset saved as drug_indications")

        result_array.append("Full dataframes are available for viewing using the command: print {dataframe_name}. See right side panel for more information.")
        result_array.append("The suffix for the drug-disease interaction pair is: %s" % query_name)
        return result_array
        
_DrugDisease = DrugDisease()


# class DrugDiseaseMulti(IrisCommand):
#     # what iris will call the command + how it will appear in a hint
#     title = "Can you find the mechanism of action of this list of drug-disease pairs"
    
#     # give an example for iris to recognize the command
#     examples = ["multiple drugs and diseases", "how list of drugs works", "how multiple drug works on disease", "how do these drugs affect these diseases", "how does this list of {drug_list} affect {disease_list}" ]

#     # type annotations for each command argument, to help Iris cosllect missing values from a user
#     argument_types = {"drug_list":t.List("What drugs do you want to analyze?"), "disease_list":t.List("What diseases do you want to analyze? List must be the same length as the list of drugs.")}
    
#     # core logic of the command
#     def command(self, drug_list, disease_list):
#         # import numpy
#         # return numpy.random.randint(100)
#         assert(len(drug_list) == len(disease_list))
#         results = []
#         for drug, disease in zip(drug_list, disease_list):
#             result = run_test.run_drug_single(drug, disease)
#             result_txt = ' '.join(('Answer Found:', result[0], 'Answer:', result[1]))
#             results.append(' '.join((drug, disease, result_txt)))
#         return '\t'.join(results)
        
#     # wrap the output of a command to display to user
#     # by default this will be an identity function
#     # each element of the list defines a separate chat bubble
#     def explanation(self, results):

#         # return ["Book a flight to DC to get the answer. Here's your lucky number", result]

#         return results.split('\t')

# _DrugDiseaseMulti = DrugDiseaseMulti()



# class DrugDiseaseCSV(IrisCommand):
#     # what iris will call the command + how it will appear in a hint
#     title = "Can you find the mechanism of action of these drug-disease pairs from csv {paired_csv} and saved to {saved_dir}"
    
#     # give an example for iris to recognize the command
#     examples = ["multiple drugs and diseases from csv", "csv of drug disease matches", "loaded drug disease match", "find the mechanism of action of drug-disease pairs from csv {paired_csv} and saved to {saved_dir}"]

#     # type annotations for each command argument, to help Iris cosllect missing values from a user
#     argument_types = {"paired_csv":t.File("What drugs do you want to analyze?"), "saved_dir":t.String("Where do you want to save the results")}

#     # core logic of the command
#     def command(self, paired_csv, saved_dir):
#         import pandas as pd
#         new_df = iris_objects.IrisDataframe(None, empty=True)
#         new_df.df = pd.read_csv(paired_csv.path, sep='\t')
#         # print('read in file')
#         self.iris.add_to_env('drug_disease_list', new_df)
#         # print('added to environment')

#         from ncats.scripts import run_test
#         print('path', paired_csv.path)
#         print('storage_dir', saved_dir)
#         run_test.run_drug_multi(paired_csv.path, storage_dir = saved_dir, iterate_until_found=False)
#         return ' '.join(("written results to:", saved_dir))
        
#     # wrap the output of a command to display to user
#     # by default this will be an identity function
#     # each element of the list defines a separate chat bubble
#     def explanation(self, results):

#         return results

# _DrugDiseaseCSV = DrugDiseaseCSV()


