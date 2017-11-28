
import os
import pandas
from subprocess import call
import argparse
from shutil import copy

# Set relative paths for data and scripts
overall_path = os.path.abspath(os.path.dirname(__file__))
pubmed_data_path = os.path.join(overall_path, 'ncats/PubMed/')
scripts_dir = os.path.join(overall_path, 'ncats/scripts/Q1/')

# This function attempts to answer the question:
# Do any genetic diseases protect against this condition
def Q1_query(condition):
    # Create a blank object to store results in
    ss_object = SemanticSimilarityResults(condition)
    try:
        # Compute semantic similarity
        results = semantic_similarity(condition)
        if results is None:
            raise Exception("Error computing semantic similarity")
        ss_object.similarities = results
    except Exception as e:
        ss_object.error = e
    return ss_object

# This function generates word clouds to study the semantic relationship between two provided terms
def investigate_similarity(disease_a, disease_b):
    # Set a blank object to store results in
    similarity = Cooccurrence(disease_a, disease_b)
    try:
        # Search for cooccurrences of the diseases in pubmed.  Set results to blank object
        cooccurrence = cooccurrence_seach(disease_a, disease_b)
        similarity.frequency_word_cloud = cooccurrence.frequency_word_cloud
        similarity.tfidf_word_cloud = cooccurrence.tfidf_word_cloud
        similarity.sentence_file = cooccurrence.sentence_file
        similarity.fetch_sentences()
    except Exception as e:
        print("Couldn't determine cooccurrence in pubmed")

    try:
        # Create commonality and comparison word clouds
        commonality, comparison = semantic_similarity_word_clouds(disease_a, disease_b)
        similarity.commonality_word_cloud = commonality
        similarity.comparison_word_cloud = comparison
    except Exception as e:
        print("Could not compute commonality or comparison word clouds")
    return similarity

# Search PubMed for abstracts that mention both diseases
# This is done using an R script.  The function is called from the shell and result are printed to a file.
# Pre-existing results are checked for and returned if they exist.
def cooccurrence_seach(disease_a, disease_b):
    # Check if the results already exist
    sentences_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_sentences.txt" % (clean_query(disease_a), clean_query(disease_b)))
    freq_word_cloud_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_frequency.png" % (
    clean_query(disease_a), clean_query(disease_b)))
    tfidf_word_cloud_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_tfidf.png" % (
    clean_query(disease_a), clean_query(disease_b)))
    if os.path.isfile(sentences_file) & os.path.isfile(freq_word_cloud_file) & os.path.isfile(tfidf_word_cloud_file):
        return Cooccurrence(disease_a, disease_b, sentences_file, freq_word_cloud_file, tfidf_word_cloud_file)
    # Run the script
    script = os.path.join(scripts_dir, "search_cooccurence.R")
    cmd = "Rscript %s -a \"%s\" -b \"%s\" --data_dir %s -p %s" % (script, disease_a, disease_b, pubmed_data_path, pubmed_data_path + "/cooccurence")
    print(cmd)
    call(cmd, shell=True)
    # Collect the results
    if os.path.isfile(sentences_file) & os.path.isfile(freq_word_cloud_file) & os.path.isfile(tfidf_word_cloud_file):
        return Cooccurrence(disease_a, disease_b, sentences_file, freq_word_cloud_file, tfidf_word_cloud_file)

# Run semantic similarity
# Compute the cosine similarity between the tf-idf vectors for the query disease against a set of genetic diseases
# This is done using an R script run from the shell
# Pre-existing results are checked for and returned if they exist
def semantic_similarity(disease):
    print("Checking semantic similarity")
    # Check whether this disease has been run before
    output_file = os.path.join(pubmed_data_path, "similarities.%s.txt" % clean_query(disease))
    if os.path.isfile(output_file):
        # Output already exists.  Return it
        print("Found existing file")
        return semantic_similarity_results(output_file)
    print("No precomputed results found.  Computing.")

    # Run the R script and check for the output files
    ss_script = os.path.join(scripts_dir, "semantic_similarity.R")
    cmd = "Rscript %s --query \"%s\" --data_dir %s" % (ss_script, disease, pubmed_data_path)
    call(cmd, shell=True)

    # If nothing there, return failure
    if os.path.isfile(output_file):
        # Output already exists.  Return it
        print("Semantic similarity complete")
        return semantic_similarity_results(output_file)

# Function to iteratively create word clouds.  Not used currently.
def batch_word_clouds(disease_a, diseases):
    commonality_clouds = []
    comparison_clouds = []
    for d in diseases:
        c1, c2 = semantic_similarity_word_clouds(disease_a, d)
        commonality_clouds.append(c1)
        comparison_clouds.append(c2)
    return(commonality_clouds, comparison_clouds)

# Generate word clouds
# Call an R function to create commonality and comparison word clouds
def semantic_similarity_word_clouds(disease_a, disease_b):
    print("Generating word clouds")
    # Check for the existence of the pngs
    commonality_file = "word_cloud.%s.%s.commonality.png" % (clean_query(disease_a), clean_query(disease_b))
    commonality_path = os.path.join(pubmed_data_path, commonality_file)

    comparison_file = "word_cloud.%s.%s.comparison.png" % (clean_query(disease_a), clean_query(disease_b))
    comparison_path = os.path.join(pubmed_data_path, comparison_file)

    if os.path.isfile(commonality_path) and os.path.isfile(comparison_path):
        return commonality_path, comparison_path

    # Check that there is an RDS file for a run of disease_a,
    rds_file = "tfidf.%s.rds" % clean_query(disease_a)
    rds_path = os.path.join(pubmed_data_path, rds_file)
    if not os.path.isfile(rds_path):
        semantic_similarity(disease_a)

    # Run the R script
    wc_script = os.path.join(scripts_dir, "comparison_clouds.R")
    prefix = pubmed_data_path + "/word_cloud"
    cmd = 'Rscript %s -f %s -a \"%s\" -b \"%s\" -p %s' % (wc_script, rds_path, disease_a, disease_b, prefix)
    print(cmd)
    call(cmd, shell=True)

    # Check for the existence of the pngs
    # Return the pngs
    if os.path.isfile(commonality_path) and os.path.isfile(comparison_path):
        return commonality_path, comparison_path
    print("Word clouds not found!")
    return None, None

# Return a "clean" version of the disease.  Replace spaces with underscores and remove commas
def clean_query(query):
    query = query.replace(" ", "_")
    query = query.replace(",", "")
    return query

# Process the results from a semantic similarity run and return a pandas data frame
def semantic_similarity_results(path):
    data = pandas.read_csv(path, sep=",")
    data.columns = ['Query', 'Genetic_Disease', 'Score']
    return data

# A results object for the semantic similarity query.  Makes error reporting to Iris possible and stores
# all file paths, etc. in a single place.
class SemanticSimilarityResults(object):
    def __init__(self, query_disease):
        self.query_disease = query_disease
        self.clean_q = clean_query(query_disease)
        # A string for any results encountered
        self.error = None
        # Pandas data frame with top 10 resutls
        self.similarities = None
        # List of paths to word clouds
        self.commonality_clouds = []
        self.comparison_clouds = []

    def top_similarities(self, n=10):
        # Return top ten most similar diseases by cosine similarity
        try:
            return self.similarities.nlargest(n, 'Score')
        except Exception as e:
            return None

    def top_disease_matches(self, n=10):
        try:
            top = self.top_similarities(n)
            if top is not None:
                return top['Genetic_Disease'].tolist()
            return None

        except Exception as e:
            return None

    def get_commonality_cloud(self, disease):
        for c in self.commonality_clouds:
            if not c:
                continue
            current = c.split("/")[-1]
            d_b = current.split(".")[2]
            if disease == d_b:
                return c
        return None

    def get_comparison_cloud(self, disease):
        for c in self.comparison_clouds:
            if not c:
                continue
            current = c.split("/")[-1]
            d_b = current.split(".")[2]
            if disease == d_b:
                return c
        return None

    def print_top_matches(self, n=10):
        top = self.top_similarities(n)
        if top is not None:
            print(top[['Genetic_Disease', 'Score']])
        else:
            print("No results found")



    # Print a summary of the results
    def print_summary(self, n=10):
        print()
        print('-------------------------------------------------')
        print('Search results')
        print('The following genetic diseases are the top candidates for your query: %s' % self.query_disease)
        print(self.print_top_matches(n))

    def write_results(self, file=None, n=10):
        top = self.top_similarities(n)
        if file is None:
            file = "%s.top_genetic_disease.csv" % self.clean_q
        if top is not None:
            top.to_csv(file)
        return file


# A cooccurrence object to store results of semantic relationship searches
class Cooccurrence(object):
    def __init__(self, a, b, sentence_file=None, frequency_word_cloud=None, tfidf_word_cloud=None):
        self.disease_a = a
        self.disease_b = b
        self.frequency_word_cloud = frequency_word_cloud
        self.tfidf_word_cloud = tfidf_word_cloud
        self.sentence_file = sentence_file
        self.sentences = []
        self.commonality_word_cloud = None
        self.comparison_word_cloud = None
        self.sentence_df = None
        self.word_cloud = None
        self.error = None

    # Retrieve sentences from a file
    def fetch_sentences(self):
        data = pandas.read_csv(self.sentence_file, sep="\t")
        self.sentence_df = data
        self.sentences = data['word'].tolist()

    # Get a subset of sentences as a data frame
    def top_sentence_df(self, n=10):
        try:
            sentences = [('Sentences', self.sentences[0:n])]
            if len(sentences) == 1:
                return None
            print(sentences)
            return pandas.DataFrame.from_items(sentences)
        except Exception as e:
            return None

    def print_sentences(self, n=10):
        self.fetch_sentences()
        if self.sentences is not None:
            print("Sample of incidences of terms cooccurring in PubMed in the same sentence:")
            for s in self.sentences[0:10]:
                print('"%s"' % s)
                print()


    # Print a summary of the results
    def print_summary(self):
        print()
        print('-------------------------------------------------')
        print('Search results for %s and %s' % (self.disease_a, self.disease_b))

        self.print_sentences()

        if self.comparison_word_cloud is not None:
            print("Comparison word cloud generated: %s" % os.path.basename(self.comparison_word_cloud))

        if self.commonality_word_cloud is not None:
            print("Comparison word cloud generated: %s" % os.path.basename(self.commonality_word_cloud))

        if self.frequency_word_cloud is not None:
            print("Comparison word cloud generated: %s" % os.path.basename(self.frequency_word_cloud))

        if self.sentence_file is not None:
            print("File containing PubMed sentences written: %s" % os.path.basename(self.sentence_file))


    def copy_files(self, path=None):
        if path is None:
            path = os.path.dirname(os.path.realpath(__file__))
        if self.comparison_word_cloud is not None:
            copy(self.comparison_word_cloud, path)

        if self.commonality_word_cloud is not None:
            copy(self.commonality_word_cloud, path)

        if self.frequency_word_cloud is not None:
            copy(self.frequency_word_cloud, path)

        if self.sentence_file is not None:
            copy(self.sentence_file, path)

def semantic_sim(args):
    results = Q1_query(args.condition)
    results.print_summary()

    file = None
    if args.output is not None:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        file = os.path.join(args.output, "%s.top_genetic_disease.csv" % results.clean_q)

    path = results.write_results(file)
    print("Results written to %s" % path)
    # move results to a new results directory


def similarity(args):
    results = investigate_similarity(args.a, args.b)
    results.print_summary()
    # todo move results to a new results directory

    if args.output is not None:
        if not os.path.exists(args.output):
            os.makedirs(args.output)

    results.copy_files(args.output)

    # todo clean up print summary

def parse_command_line():
    parser = argparse.ArgumentParser(
        description='NCATS Question 1 Command Line Interface')

    subparsers = parser.add_subparsers()

    semantic_parser = subparsers.add_parser('protective', help="Search for genetic diseases that might be protective of the query condition")
    semantic_parser.add_argument("condition", help="Condition for which to find protective genetic diseases")
    semantic_parser.add_argument("--output", default=None, help="Set a results directory ")
    semantic_parser.add_argument("--debug", default=False, action='store_true', help="Print debugging messages")
    semantic_parser.set_defaults(func=semantic_sim)

    similarity_parser = subparsers.add_parser('relationship', help="Search for the relationship between two diseases.  "
                                                                   "Ideal for understanding relationships returned by "
                                                                   "the protective function.")
    similarity_parser.add_argument("a", help="First disease. Use quotes if the disease name has more than one word.")
    similarity_parser.add_argument("b", help="Second disease. Use quotes if the disease name has more than one word.")
    similarity_parser.add_argument("--output", default=None, help="Set a results directory ")
    similarity_parser.add_argument("--debug", default=False, action='store_true', help="Print debugging messages")
    similarity_parser.set_defaults(func=similarity)

    options = parser.parse_args()
    if not hasattr(options, 'func'):
        parser.print_help()
        exit(0)

    options.func(options)

if __name__ == "__main__":
    options = parse_command_line()

