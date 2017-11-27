from iris import state_types as t
from iris import IrisCommand
from app.user_functions.Q1_main import investigate_similarity
from iris import iris_objects
import os
import pandas


class InvestigateSimilarity(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "How are {disease_a} and {disease_b} semantically related?"
    # give an example for iris to recognize the command
    examples = ["How are {disease_a} and {disease_b} semantically related?",
                "What is the relationship between {disease_a} and {disease_b}?",
                "Do {disease_a} and {disease_b} appear together in PubMed?"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"disease_a": t.String("What is the first disease?"),
                      "disease_b": t.String("What is the second disease?")}

    #     # type annotations for each command argument, to help Iris collect missing values from a user
    #     argument_types = {"drug":t.String("Okay, a couple more questions to set up this task. For confirmation: What is the drug you want to analyze?"),
    #                         "disease":t.String("What is the disease you want to analyze?"),
    #                         "bool_image":t.YesNo("Would you like to visual the results as a diagram?",
    #                                     yes=True, no=False),
    #                         "bool_other_disease":t.YesNo("Would you like to know other diseases that can be treated by this drug?",
    #                                     yes=True, no=False),
    #                         "bool_pubmed":t.YesNo("Would you like to get the list of pubmed IDs for reference?",
    #                                     yes=True, no=False)
    #
    #                         }

    # ,"genetic_disease":t.String("What is the genetic disease do you think it might link to? If unknown, type none")}

    # core logic of the command
    def command(self, disease_a, disease_b):
        # Run the query
        results = investigate_similarity(disease_a, disease_b)

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

        if result is None:
            result_array.append('There was an error processing your request')
            return result_array


        if result.error is not None:
            result_array.append('There was an error processing your request')
            return result_array

        sentence_df = result.top_sentence_df()
        if sentence_df is None:
            result_array.append('The two conditions entered do not appear together in the same sentence in PubMed')
        else:
            result_array.append('Here are some examples of the two diseases appearing in the same sentence in PubMed')
            result_array.append(sentence_df)
            sentence_df = iris_objects.IrisDataframe(data=result.top_sentence_df())
            self.iris.add_to_env('sentences', sentence_df)



        if result.commonality_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a commonality word cloud. This shows terms enriched between the two diseases '
                                'in PubMed')
            os.system("open " + result.commonality_word_cloud)

        if result.comparison_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a comparison word cloud. This shows terms enriched in each disease separately '
                                'in PubMed.')
            os.system("open " + result.comparison_word_cloud)

        if result.frequency_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a cooccurrence word cloud. This shows terms enriched in PubMed where the two '
                                'conditions appear together in the same abstract.')
            os.system("open " + result.frequency_word_cloud)

        if len(result_array) < 1:
            result_array.append("No similarity metrics available")

        return result_array


_InvestigateSimilarity = InvestigateSimilarity()