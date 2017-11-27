from iris import state_types as t
from iris import IrisCommand
from app.user_functions.Q1_main import Q1_query
from iris import iris_objects
import os

class GeneticConditionDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "What genetic disease protects against {condition}?"
    # give an example for iris to recognize the command
    examples = ["what protects against {condition}", "what genetic conditions might offer protection against {condition} and why", "protective mechanism of condition against disease", "protective condition", "protection against disease"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"condition":t.String("Just for confirmation: What is the condition do you want to analyze?")} 
                        # ,"genetic_disease":t.String("What is the genetic disease do you think it might link to? If unknown, type none")}
    
    # core logic of the command
    def command(self, condition):
        # Run the query
        results = Q1_query(condition)

        return results

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        """"# results is an object with       # Pandas data frame with top 10 resutls
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        """
        result_array = []
        if result.error is not None:
            result_array.append('There was an error processing your request')
            return result_array

        result_array.append('The following genetic diseases are the most semantically similar to your query')
        # adds the table to results
        similarities_df = iris_objects.IrisDataframe(data=result.top_similarities())
        self.iris.add_to_env('similarities', similarities_df)
        result_array.append(similarities_df)

        # display image (first one)
        #if len(result.commonality_clouds > 0):
        #    os.system("open " + result.commonality_clouds[0])


        return result_array


_GeneticConditionDisease = GeneticConditionDisease()