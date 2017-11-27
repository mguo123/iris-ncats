
import os
import pandas
from subprocess import call


overall_path = os.path.abspath(os.path.dirname(__file__))
pubmed_data_path = os.path.join(overall_path, 'ncats/PubMed/')
scripts_dir = os.path.join(overall_path, 'ncats/scripts/Q1/')

# This function attempts to answer the question:
# Do any genetic diseases protect against this condition
def Q1_query(condition):
    ss_object = SemanticSimilarityResults(condition)
    try:
        # Compute semantic similarity
        print("Computing similarity")
        results = semantic_similarity(condition)
        if results is None:
            raise Exception("Error computing semantic similarity")
        ss_object.similarities = results

        # Get word cloud paths
        #commonality, comparison = batch_word_clouds(condition, ss_object.top_disease_matches(3))
        #ss_object.commonality_clouds = commonality
        #ss_object.comparison_clouds = comparison


    except Exception as e:
        ss_object.error = e

    #if ss_object.error is None:
        #print("FINAL RESULTS")
        #print(ss_object.top_similarities())
        #print(ss_object.commonality_clouds)
        #print(ss_object.comparison_clouds)
    #else:
    #    print(ss_object.error)
    return ss_object

def investigate_similarity(disease_a, disease_b, word_cloud=None):
    # 1. Cooccurrence search
    similarity = cooccurence_seach(disease_a, disease_b)

    # 2. Commonality word cloud

    if word_cloud == "commonality":
        print("Creating commonality wordcloud")
        commonality, comparison = semantic_similarity_word_clouds(disease_a, disease_b)
        similarity.word_cloud = commonality
    elif word_cloud == "comparison":
        print("Creating comparison wordcloud")
        commonality, comparison = semantic_similarity_word_clouds(disease_a, disease_b)
        similarity.word_cloud = comparison
    elif word_cloud == "cooccurrence":
        print("Creating cooccurrence wordcloud")
        similarity.word_cloud = similarity.frequency_word_cloud
    else:
        similarity.word_cloud = None

    return similarity


def cooccurence_seach(disease_a, disease_b):
    # Check if the results already exist
    sentences_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_sentences.txt" % (clean_query(disease_a), clean_query(disease_b)))
    freq_word_cloud_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_frequency.png" % (
    clean_query(disease_a), clean_query(disease_b)))
    tfidf_word_cloud_file = os.path.join(pubmed_data_path, "cooccurence.%s.%s.cooccurence_tfidf.png" % (
    clean_query(disease_a), clean_query(disease_b)))

    if os.path.isfile(sentences_file) & os.path.isfile(freq_word_cloud_file) & os.path.isfile(tfidf_word_cloud_file):

        return Cooccurence(disease_a, disease_b, sentences_file, freq_word_cloud_file, tfidf_word_cloud_file)

    # Run the script
    script = os.path.join(scripts_dir, "search_cooccurence.R")
    cmd = "Rscript %s -a \"%s\" -b \"%s\" --data_dir %s -p %s" % (script, disease_a, disease_b, pubmed_data_path, pubmed_data_path + "/cooccurence")
    print(cmd)
    call(cmd, shell=True)

    # Collect the results
    if os.path.isfile(sentences_file) & os.path.isfile(freq_word_cloud_file) & os.path.isfile(tfidf_word_cloud_file):
        return Cooccurence(disease_a, disease_b, sentences_file, freq_word_cloud_file, tfidf_word_cloud_file)
    print(sentences_file)
    print(freq_word_cloud_file)
    print(tfidf_word_cloud_file)


# Run semantic similarity
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


def batch_word_clouds(disease_a, diseases):
    print(disease_a)
    print(diseases)

    commonality_clouds = []
    comparison_clouds = []

    for d in diseases:
        c1, c2 = semantic_similarity_word_clouds(disease_a, d)
        commonality_clouds.append(c1)
        comparison_clouds.append(c2)
    return(commonality_clouds, comparison_clouds)

    # Check if word clouds already exists

    # If not, generate the word clouds



# Generate word clouds
def semantic_similarity_word_clouds(disease_a, disease_b):
    print("Generating word clouds")
    # Check for the existence of the pngs
    commonality_file = "word_cloud.%s.%s.commonality.png" % (clean_query(disease_a), clean_query(disease_b))
    commonality_path = os.path.join(pubmed_data_path, commonality_file)

    comparison_file = "word_cloud.%s.%s.comparison.png" % (clean_query(disease_a), clean_query(disease_b))
    comparison_path = os.path.join(pubmed_data_path, comparison_file)

    print(commonality_path)
    #print(comparison_path)


    if os.path.isfile(commonality_path) and os.path.isfile(comparison_path):
        return commonality_path, comparison_path

    # Check that there is an RDS file for a run of disease_a,
    rds_file = "tfidf.%s.rds" % clean_query(disease_a)
    rds_path = os.path.join(pubmed_data_path, rds_file)
    if not os.path.isfile(rds_path):
        semantic_similarity(disease_a)

    # Bonus points for checking whether disease b is a genetic disease
    # maybe later

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
    print("not found!")
    return None, None

def clean_query(query):
    query = query.replace(" ", "_")
    query = query.replace(",", "")
    return query

# Process the results from a semantic similarity run and return a pandas data frame
def semantic_similarity_results(path):
    data = pandas.read_csv(path, sep=",")
    return data

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
            return self.similarities.nlargest(n, 'cos_sim').head(n=n)
        except Exception as e:
            return None

    def top_disease_matches(self, n=10):
        try:
            top = self.top_similarities(n)
            if top is not None:
                return top['disease_b'].tolist()
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

    def print_summary(self):
        print(self.top_disease_matches())
        print(self.commonality_clouds)

class Cooccurence(object):
    def __init__(self, a, b, sentence_file, frequency_word_cloud, tfidf_word_cloud):
        self.disease_a = a
        self.disease_b = b
        self.frequency_word_cloud = frequency_word_cloud
        self.tfidf_word_cloud = tfidf_word_cloud
        self.sentence_file = sentence_file
        self.sentences = []
        self.fetch_sentences()
        self.commonality_word_cloud = None
        self.comparison_word_cloud = None
        self.sentence_df = None
        self.word_cloud = None
        self.error = None

    def fetch_sentences(self):
        data = pandas.read_csv(self.sentence_file, sep="\t")
        self.sentence_df = data
        self.sentences = data['word'].tolist()

    def top_sentence_df(self, n=10):
        try:
            sentences = [('Sentences', self.sentences[0:n])]
            print(sentences)
            return pandas.DataFrame.from_items(sentences)
        except Exception as e:
            return None

    def print_summary(self):
        print(self.word_cloud)
        print(self.top_sentence_df())

if __name__ == "__main__":
    condition="osteoporosis"
    Q1_query(condition)

    #condition="rabies"
    #Q1_query(condition)
    c = investigate_similarity("rubella", "grant syndrome", word_cloud="commonality")

    c.print_summary()


    #for s in c.sentences[1:10]:
    #    print(s)