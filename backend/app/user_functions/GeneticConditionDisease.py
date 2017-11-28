from iris import state_types as t
from iris import IrisCommand
from app.user_functions.Q1_main import Q1_query
from iris import iris_objects
import os

class GeneticConditionDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "What genetic disease might protect against {condition}?"
    # give an example for iris to recognize the command
    examples = ["what protects against {condition}",
                "what genetic conditions might offer protection against {condition} and why",
                "protective mechanism of condition against disease",
                "protective condition",
                "protection against disease",
                "do any genetic diseases protect against {condition}"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"condition":t.String("What is the condition do you want to analyze?")}

    # core logic of the command
    def command(self, condition):
        # Run the query
        results = Q1_query(condition)

        return [condition, results]

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        """"
        results is an object with Pandas data frame with top 10 resutls
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        """
        condition, results = result
        
        # make the df_name:
        df_name = 'similarities_' + condition[:min(len(condition), 5)]
        # remove spaces and make lowercase
        df_name = df_name.replace(" ", "")
        df_name = df_name.lower()

        result_array = []
        if results.error is not None:
            result_array.append('There was an error processing your request for %s' % condition)
            return result_array

        # adds the table to results
        similarities = results.top_similarities()
        if similarities is None:
            result_array.append('No similarities could be computed')
        else:

            similarities_df = iris_objects.IrisDataframe(data=similarities)
            # adds the table to results
            self.iris.add_to_env(df_name, similarities_df)
            result_array.append(similarities_df)

        # display image (first one)
        #if len(results.commonality_clouds > 0):
        #    os.system("open " + results.commonality_clouds[0])


        return result_array


_GeneticConditionDisease = GeneticConditionDisease()