from iris import state_types as t
from iris import IrisCommand
from app.user_functions.Q1_main import investigate_similarity
from iris import iris_objects
import os

class InvestigateSimilarity(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "How are disease_a and disease_b semantically related? "
    # give an example for iris to recognize the command
    examples = ["semantically",
                "What is the relationship between {disease_a} and {disease_b}?", "semantic relationship", "semantic relationship of diseases",
                "Do {disease_a} and {disease_b} appear together in PubMed"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"disease_a": t.String("What is the first disease?"),
                      "disease_b": t.String("What is the second disease?")}

    # core logic of the command
    def command(self, disease_a, disease_b):
        # Run the query
        results = investigate_similarity(disease_a, disease_b)

        return [results, disease_a, disease_b]

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        """"
        results -  an object with Pandas data frame with top 10 resutls that is called with top_sentences_df
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        """
        [results, disease_a, disease_b] = result
        # create name df
        df_name = 'sentences_' + disease_a[:min(len(disease_a), 3)] + '_' + disease_b[:min(len(disease_b), 3)]
        df_name = df_name.replace(" ", "")
        df_name = df_name.lower()


        result_array = []

        if results is None:
            result_array.append('There was an error processing your request')
            return result_array

        if results.error is not None:
            result_array.append('There was an error processing your request')
            return result_array

        sentence_df = results.top_sentence_df()

        if sentence_df is None:
            result_array.append('The two conditions entered do not appear together in the same sentence in PubMed')
        else:
            result_array.append('Here are some examples of the two diseases appearing in the same sentence in PubMed')
            result_array.append(sentence_df)
            sentence_df = iris_objects.IrisDataframe(data=results.top_sentence_df())
            self.iris.add_to_env(df_name, sentence_df)


        if results.commonality_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a commonality word cloud. This shows terms enriched between the two diseases '
                                'in PubMed. Stored at %s' % results.commonality_word_cloud)
            os.system("open " + results.commonality_word_cloud)

        if results.comparison_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a comparison word cloud. This shows terms enriched in each disease separately '
                                'in PubMed. Stored at %s' % results.comparison_word_cloud)
            os.system("open " + results.comparison_word_cloud)

        if results.frequency_word_cloud is not None:
            # display image (first one)
            result_array.append('Generated a cooccurrence word cloud. This shows terms enriched in PubMed where the two '
                                'conditions appear together in the same abstract. Stored at %s' % results.frequency_word_cloud)
            os.system("open " + results.frequency_word_cloud)

        if len(result_array) < 1:
            result_array.append("No similarity metrics available")


        if len(result_array) < 1:
            result_array.append("No similarity metrics available")

        return result_array

_InvestigateSimilarity = InvestigateSimilarity()
