from iris import state_types as t
from iris import IrisCommand
from iris import iris_objects

from app.user_functions.ncats.tools.GNBR_api import GNBR_api 


# from iris import iris_objects
import os
import numpy as np


class QueryDiseaseTreatments(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "What drugs treat {disease}?"
    # give an example for iris to recognize the command

    examples = ["What are the treatments for {disease}?",
                "How is {disease} treated?"]

    # type annotations for each command argument, to help Iris collect missing values from a user

    argument_types = {"disease": t.String("What disease do you want to get treatments for?")}

    # core logic of the command
    def command(self, disease):
        # Run the query
        #
        results = GNBR_api.query_treatments_for_disease(disease) # list of lists

        return [results, disease]


    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        """"# results is an object with       # Pandas data frame with top 10 resutls
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        """
        [results, disease] = result
        # make the df_name:
        df_name = 'treatments' + disease[:min(len(disease), 5)]
        # remove spaces and make lowercase
        df_name = df_name.replace(" ", "")
        df_name = df_name.lower()


        result_array = []
        if len(results) == 0:
            result_array.append('No treatments identified in GNBR')
            return result_array

        treatment_data = iris_objects.IrisDataframe(column_names=[ "Phenotype", "Mesh ID", "Frequency of Annotation"], column_types=["Text", "Text", "Text"], data=results)

        result_array.append('The following treatments were identified in GNBR for %s' % disease)
        result_array.append(treatment_data)
        return result_array


_QueryDiseaseTreatments = QueryDiseaseTreatments()


