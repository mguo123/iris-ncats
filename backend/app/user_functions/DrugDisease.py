from iris import state_types as t
from iris import IrisCommand


from iris import state_machine as sm
from iris import util as util
from iris import iris_objects

from app.user_functions.ncats.Q2.PathFx_api import run_main
from app.user_functions.ncats.Q2.Q2_main import Q2_query
import numpy as np
import os
import datetime

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
                        "bool_other_disease":t.YesNo("Would you like to know other diseases that can be affected by this drug?",
                                    yes=True, no=False)
                            
                        }
    
    # core logic of the command
    def command(self, drug, disease, bool_image, bool_pubmed, bool_other_disease):
        print('BEFORE QUERY!!!')
        # create options structure 

        options = Options(gen_image=bool_image, gen_pubmed=bool_pubmed)


        answer = Q2_query(drug, disease, options)
        print('ran Q2 query!!!')

        # Error handling in disease or drug is not found
        if isinstance(answer, str):
            answer_str = answer
            answer = {}
            answer['error'] = answer_str

        # if want to find other indications
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
        
        if "error" in result:
            return result['error']

        query_name = result['drug'][:min(len(result['drug']), 3)] + "_" + result['disease'][:min(len(result['disease']), 3)]
        query_name = '_'.join(query_name.split(' '))
        query_name = "_" + query_name.lower()
        print(query_name, "== query name")

        # Print out genes associated with drug
        query_statement = 'How does ' + result['drug'] + '(' + result['drug_id'] + ') treat ' + result['disease'] + '(' + result['disease_id'] + ').'
        result_array = ['Here are your results for: %s' % query_statement]
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
        if 'tissue_df_dis' in result:
            result_array.append('Tissues that show the most differential gene expression in the disease state are shown. Full dataset saved as tissues_{drug_disease} ')
            tissue_object_dis = iris_objects.IrisDataframe(data=result['tissue_df_dis'])
            tissue_object_dis_short = iris_objects.IrisDataframe(data=result['tissue_df_dis_short'])
            self.iris.add_to_env('tissues_disease' + query_name, tissue_object_dis)
            result_array.append(tissue_object_dis_short)
        else:
            result_array.append('No differential tissue expression in disease state detected.')
        #  # get tissue = drug
        # result_array.append('Top drug-relevant tissues in which this drug-disease interaction takes place are shown. Full dataset saved as tissues_{drug_disease} ')
        # tissue_object_drug = iris_objects.IrisDataframe(data=result['tissue_df_drug'])
        # tissue_object_drug_short = iris_objects.IrisDataframe(data=result['tissue_df_drug_short'])
        # self.iris.add_to_env('tissues_drug' + query_name, tissue_object_drug)
        # result_array.append(tissue_object_drug_short)


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
            if len(ph_genes_arr) >=4:
                ph_genes_array_all_iris = iris_objects.IrisDataframe(column_names=[ "Phenotype", "probability", "Benjamin Hochberg significance cutoff","list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_all)
                self.iris.add_to_env('drug_indications' + query_name, ph_genes_array_all_iris)
                ph_genes_array_short = [ph_genes_arr[x:x+4] for x in range(0, min(5*4,len(ph_genes_arr)),4)]
                ph_genes_array_short_iris = iris_objects.IrisDataframe(column_names=["Phenotype", "Probability", "Benjamin Hochberg significance cutoff", "list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_short)
                result_array.append(ph_genes_array_short_iris)  
            # result_array.append("Full dataset saved as drug_indications")


        # display image
        if "image_file" in result:
            result_array.append('Diagram stored in: %s' % result["image_file"])
            os.system("open " + result["image_file"])

            
        result_array.append("Full dataframes are available for viewing using the command: print {dataframe_name}. See right side panel for more information.")
        result_array.append("The suffix for the drug-disease interaction pair is: %s" % query_name)
        result_array.append("Results are also stored in: %s" % results_dir)

        return result_array
        
_DrugDisease = DrugDisease()





class DrugDiseaseMulti(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "Can you find the mechanism of action of this list of drug-disease pairs"

    # give an example for iris to recognize the command
    examples = ["multiple drugs and diseases", "how list of drugs works", "how multiple drug works on disease", "how do these drugs affect these diseases", "how does this list of {drug_list} affect {disease_list}" ]

    # type annotations for each command argument, to help Iris cosllect missing values from a user
    argument_types = {"drug_list":t.List("What drugs do you want to analyze? Please enter in drugs separated by commas."), 
                        "disease_list":t.List("What diseases do you want to analyze? Every pairwise combination of drugs and diseases will be computed."),
                        "bool_image":t.YesNo("Would you like to save visual representations of the drug-disease combination?",
                                    yes=True, no=False),
                        "bool_pubmed":t.YesNo("Would you like to get the list of pubmed IDs for reference for each query?",
                                    yes=True, no=False),
                        "bool_other_disease":t.YesNo("Would you like to know other diseases that can be affected by the given drug?",
                                    yes=True, no=False)
                            
                        }

    # core logic of the command
    def command(self, drug_list, disease_list, bool_image, bool_pubmed, bool_other_disease):
        # generate options object
        options = Options(gen_image=bool_image, gen_pubmed=bool_pubmed)

        # store the answers to each drug-disease combination as a list of dictionaries
        answer_arr = []

        for drug in drug_list:
            for disease in disease_list:

                answer = Q2_query(drug, disease, options)
                if isinstance(answer, str):
                    answer_str = answer
                    answer = {}
                    answer['error'] = answer_str

                if bool_other_disease:
                    answer["other_disease"] = run_main.find_drug_indications(drug)

                answer['drug'] = drug
                answer['disease'] = disease

                answer_arr.append(answer)
    
        return answer_arr



    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, results):

        explanation_array = ['Result tables for each query are stored in the right side bar as variables. You can view a table using the command: print {dataframe_name}_{suffix}.']
        explanation_array.append('Diagrams (if requested) and other results can be found in the results directory: %s' % results_dir)
        explanation_array.append('Suffix and variable information is displayed below')
        # iterate through every drug-disease pair
        drug_arr = []
        disease_arr = []
        worked_arr = []
        suffix_arr = []
        assoc_variables = []

        for i, result in enumerate(results):
            print('results', result['drug'], result['disease'])
            drug_arr.append(result['drug'])
            disease_arr.append(result['disease'])

            if 'error' in result:
                worked_arr.append(result['error'])
                suffix_arr.append('')
                assoc_variables.append('')
            else:
                worked_arr.append('SUCCESS')

                # get suffix information
                query_name = result['drug'][:min(len(result['drug']), 3)] + "_" + result['disease'][:min(len(result['disease']), 3)] 
                query_name = '_'.join(query_name.split(' '))
                query_name = "_" + query_name.lower() + "_" + str(i)
                suffix_arr.append(query_name)

                # get associated drug genes
                drug_gene_term_object = iris_objects.IrisDataframe(data=result['drug_genes'])
                self.iris.add_to_env('drug_genes' + query_name, drug_gene_term_object) 

        
                # get genes associated with disease
                disease_gene_term_object = iris_objects.IrisDataframe(data=result['disease_genes'])
                self.iris.add_to_env('disease_genes' + query_name, disease_gene_term_object)

                # get out signficant go terms
                go_term_object = iris_objects.IrisDataframe(data=result['GOENRICH'])
                self.iris.add_to_env('go_terms' + query_name, go_term_object)

                variable_info = ['drug_genes' + query_name, 'disease_genes' + query_name, 'go_terms' + query_name]

                # get tissue = disease
                if 'tissue_df_dis' in result:
                    variable_info.append('tissues_disease' + query_name)
                    tissue_object_dis = iris_objects.IrisDataframe(data=result['tissue_df_dis'])
                    self.iris.add_to_env('tissues_disease' + query_name, tissue_object_dis)
 
                # # get tissue = drug
                # tissue_object_drug = iris_objects.IrisDataframe(data=result['tissue_df_drug'])
                # self.iris.add_to_env('tissues_drug' + query_name, tissue_object_drug)


                if "pubmed" in result:
                    if not isinstance(result["pubmed"], str):
                        variable_info.append('pmid' + query_name)       
                        pmid_df = iris_objects.IrisDataframe(data=result["pubmed"])
                        self.iris.add_to_env('pmid' + query_name, pmid_df)
                       

                # get other possible disease 
                if "other_disease" in result:
                    ph_genes_str, drug = result["other_disease"]
                    ph_genes_arr = ph_genes_str.split('\t') # prb, BH, ph, sig_genes 
                    if len(ph_genes_arr) >=4:
                        ph_genes_array_all = [ph_genes_arr[x:x+4] for x in range(0, len(ph_genes_arr),4)]
                        ph_genes_array_all_iris = iris_objects.IrisDataframe(column_names=[ "Phenotype", "probability", "Benjamin Hochberg significance cutoff","list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_all)
                        self.iris.add_to_env('drug_indications' + query_name, ph_genes_array_all_iris)
                        variable_info.append('drug_indications' + query_name)  

                assoc_variables.append(', '.join(variable_info))   


        # Save info as an iris dataframe
        info_data = [list(x) for x in zip(drug_arr, disease_arr, worked_arr, suffix_arr, assoc_variables)]

        info_df = iris_objects.IrisDataframe(column_names=[ "Drug", "Disease", "Query Status","Suffix", "Associated Variables"], column_types=["Text", "Text", "Text", "Text", "Text"], data=info_data)
        explanation_array.append(info_df)

        return explanation_array

_DrugDiseaseMulti = DrugDiseaseMulti()



class DrugDiseaseCSV(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "Can you find the mechanism of action of these drug-disease pairs from {file}"

    # give an example for iris to recognize the command
    examples = ["multiple drugs and diseases from csv", "csv of drug disease matches", "loaded drug disease match", "find the mechanism of action of drug-disease pairs from {file} "]

    # type annotations for each command argument, to help Iris cosllect missing values from a user
    argument_types = {"file":t.File("Enter your tab-delimited file containing drug-disease pairs of connections do you want to analyze. Drugs should be in the first column and conditions in the second column"), 
                        "bool_image":t.YesNo("Would you like to save visual representations of the drug-disease combination?",
                                    yes=True, no=False),
                        "bool_pubmed":t.YesNo("Would you like to get the list of pubmed IDs for reference for each query?",
                                    yes=True, no=False),
                        "bool_other_disease":t.YesNo("Would you like to know other diseases that can be affected by the given drug?",
                                    yes=True, no=False),
                        "bool_display":t.YesNo("Would you like to display results within the Iris GUI. Please only do so if you have < 10 comparisons.",
                                    yes=True, no=False)                            
                            }

    # core logic of the command
    def command(self, file, bool_image, bool_pubmed, bool_other_disease, bool_display):
        import pandas as pd
        panda_df = pd.read_csv(file.path, sep='\t')
        iris_df = iris_objects.IrisDataframe(data=panda_df)
        
        # print('read in file')
        self.iris.add_to_env('drug_disease_list', iris_df)


        # generate options object
        task_dir = os.path.join(results_dir, datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        options = Options(gen_image=bool_image, gen_pubmed=bool_pubmed, outPath=task_dir)

        # get list of drugs and conditions
        drug_list = list(panda_df.ix[:, 0])
        disease_list = list(panda_df.ix[:, 1])

        # store within the directory as a list
        answer_arr = [task_dir]
        for drug, disease in zip(drug_list, disease_list):

            answer = Q2_query(drug, disease, options)
            if isinstance(answer, str):
                answer_str = answer
                answer = {}
                answer['error'] = answer_str

            if bool_other_disease:
                answer["other_disease"] = run_main.find_drug_indications(drug)

            answer['drug'] = drug
            answer['disease'] = disease
            answer_arr.append(answer)

        if bool_display:
            return answer_arr
        else:
            return ' '.join(("written results to:", task_dir))

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble


    #### TO DO MAKE THIS A SUBCLASS OF DRUGDISEASEMULTI####
    def explanation(self, results):
        if isinstance(results, str):
            return results

        else:
            task_dir = results.pop(0)

            explanation_array = ['Result tables for each query are stored in the right side bar as variables. You can view a table using the command: print {dataframe_name}_{suffix}.']
            explanation_array.append('Diagrams (if requested) and other results can be found in the results directory: %s' % task_dir)
            explanation_array.append('Suffix and variable information is displayed below')
            # iterate through every drug-disease pair
            drug_arr = []
            disease_arr = []
            worked_arr = []
            suffix_arr = []
            assoc_variables = []

            for result in results:
                print('results', result['drug'], result['disease'])
                drug_arr.append(result['drug'])
                disease_arr.append(result['disease'])

                if 'error' in result:
                    worked_arr.append(result['error'])
                    suffix_arr.append('')
                    assoc_variables.append('')
                else:
                    worked_arr.append('SUCCESS')

                    # get suffix information
                    query_name = result['drug'][:min(len(result['drug']), 3)] + "_" + result['disease'][:min(len(result['disease']), 3)]
                    query_name = '_'.join(query_name.split(' '))
                    query_name = "_" + query_name.lower()
                    suffix_arr.append(query_name)

                    # get associated drug genes
                    drug_gene_term_object = iris_objects.IrisDataframe(data=result['drug_genes'])
                    self.iris.add_to_env('drug_genes' + query_name, drug_gene_term_object) 

            
                    # get genes associated with disease
                    disease_gene_term_object = iris_objects.IrisDataframe(data=result['disease_genes'])
                    self.iris.add_to_env('disease_genes' + query_name, disease_gene_term_object)

                    # get out signficant go terms
                    go_term_object = iris_objects.IrisDataframe(data=result['GOENRICH'])
                    self.iris.add_to_env('go_terms' + query_name, go_term_object)


                    variable_info = ['drug_genes' + query_name, 'disease_genes' + query_name, 'go_terms' + query_name]

                    # get tissue = disease
                    if 'tissue_df_dis' in result:
                        variable_info.append('tissues_disease' + query_name)
                        tissue_object_dis = iris_objects.IrisDataframe(data=result['tissue_df_dis'])
                        self.iris.add_to_env('tissues_disease' + query_name, tissue_object_dis)
     

                    if "pubmed" in result:
                        if not isinstance(result["pubmed"], str):
                            variable_info.append('pmid' + query_name)       
                            pmid_df = iris_objects.IrisDataframe(data=result["pubmed"])
                            self.iris.add_to_env('pmid' + query_name, pmid_df)
                           

                    # get other possible disease 
                    if "other_disease" in result:
                        ph_genes_str, drug = result["other_disease"]
                        ph_genes_arr = ph_genes_str.split('\t') # prb, BH, ph, sig_genes 
                        if len(ph_genes_arr) >=4:
                            ph_genes_array_all = [ph_genes_arr[x:x+4] for x in range(0, len(ph_genes_arr),4)]
                            ph_genes_array_all_iris = iris_objects.IrisDataframe(column_names=[ "Phenotype", "probability", "Benjamin Hochberg significance cutoff","list of genes"], column_types=["Text", "Text", "Text", "Text"], data=ph_genes_array_all)
                            self.iris.add_to_env('drug_indications' + query_name, ph_genes_array_all_iris)
                            variable_info.append('drug_indications' + query_name)  

                    assoc_variables.append(', '.join(variable_info))   

            # Save info as an iris dataframe
            info_data = [list(x) for x in zip(drug_arr, disease_arr, worked_arr, suffix_arr, assoc_variables)]
            info_df = iris_objects.IrisDataframe(column_names=[ "Drug", "Disease", "Query Status","Suffix", "Associated Variables"], column_types=["Text", "Text", "Text", "Text", "Text"], data=info_data)
            explanation_array.append(info_df)

            return explanation_array

_DrugDiseaseCSV = DrugDiseaseCSV()


