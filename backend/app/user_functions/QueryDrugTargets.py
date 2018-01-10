from iris import state_types as t
from iris import IrisCommand
from app.user_functions.ncats.scripts.Q1_main import Q1_query
from app.user_functions.ncats.tools.orange_apis import QueryUniprot
from app.user_functions.ncats.tools.orange_apis import QueryChEMBL


# from iris import iris_objects
import os
import numpy as np


class QueryDrugTargets(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "What genes does drug target?"
    # give an example for iris to recognize the command

    examples = ["What are the targets of {drug}",
                "What are the drug targets of {drug}",
                "Search ChEMBL for targets of {drug}"]

    # type annotations for each command argument, to help Iris collect missing values from a user

    argument_types = {"drug": t.String("What drug do you want to get gene targets for?")}

    # core logic of the command
    def command(self, drug):
        # Run the query
        #
        results = self.fetch_gene_targets(drug)

        return [results, drug]


    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        """"# results is an object with       # Pandas data frame with top 10 resutls
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        """
        [results, drug] = result
        # make the df_name:
        df_name = 'target_genes_' + drug[:min(len(drug), 5)]
        # remove spaces and make lowercase
        df_name = df_name.replace(" ", "")
        df_name = df_name.lower()


        result_array = []
        if len(results) == 0:
            result_array.append('No drug targets identified in ChEMBL')
            return result_array

        gene_results_str = ', '.join(results)

        result_array.append('The following genes were identified as drug targets in ChEMBL for %s' % drug)
        result_array.append(gene_results_str)
        return result_array

    def fetch_gene_targets(self,drug):
        # Fetch the uniprot ids for the drug targets from Chembl
        try:
            print("Fetching uniprot ids from ChEMBL for drug targets")
            uniprot_ids = QueryChEMBL.QueryChEMBL.get_target_uniprot_ids_for_drug(drug)

            print("Fetching gene names for uniprot ids")
            gene_results = np.array([])
            count = 0
            for u in uniprot_ids:
                gene = QueryUniprot.QueryUniprot.uniprot_id_to_gene_name(u)
                gene_results = np.append(gene_results,np.array(list(gene)))
                if count >= 5:
                    remaining = len(uniprot_ids) - count
                    if remaining < 1:
                        break
                    gene_results = np.append(gene_results,np.array(list("and %s more" % remaining)))
                    break
                count += 1

            print(gene_results)
            return gene_results
        except Exception as e:
            print("Error fetching gene targets: %s" % e)
            return 'Error'


_QueryDrugTargets = QueryDrugTargets()

# if __name__ == "__main__":
#     # q = QueryDrugTargets()
#     print(fetch_gene_targets('clothiapine'))
