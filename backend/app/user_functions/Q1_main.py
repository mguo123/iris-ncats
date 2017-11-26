
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
        commonality, comparison = batch_word_clouds(condition, ss_object.top_disease_matches())
        ss_object.commonality_clouds = commonality
        ss_object.comparison_clouds = comparison


    except Exception as e:
        ss_object.error = e

    if ss_object.error is None:
        print("FINAL RESULTS")
        print(ss_object.top_similarities())
        print(ss_object.commonality_clouds)
        print(ss_object.comparison_clouds)
    else:
        print(ss_object.error)
    return ss_object


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
    commonality_file = "word_cloud.%s.%s.commonality.png" % (disease_a, disease_b)
    commonality_path = os.path.join(pubmed_data_path, commonality_file)

    comparison_file = "word_cloud.%s.%s.comparison.png" % (disease_a, disease_b)
    comparison_path = os.path.join(pubmed_data_path, comparison_file)

    if os.path.isfile(commonality_path) and os.path.isfile(comparison_path):
        return commonality_path, comparison_path

    # Check that there is an RDS file for a run of disease_a,
    rds_file = "tfidf.%s.rds" % disease_a
    rds_path = os.path.join(pubmed_data_path, rds_file)
    if not os.path.isfile(rds_path):
        semantic_similarity(disease_a)

    # Bonus points for checking whether disease b is a genetic disease
    # maybe later

    # Run the R script
    wc_script = os.path.join(scripts_dir, "comparison_clouds.R")
    prefix = pubmed_data_path + "/word_cloud"
    cmd = 'Rscript %s -f %s -a %s -b %s -p %s' % (wc_script, rds_path, disease_a, disease_b, prefix)
    call(cmd, shell=True)

    # Check for the existence of the pngs
    # Return the pngs
    if os.path.isfile(commonality_path) and os.path.isfile(comparison_path):
        return commonality_path, comparison_path

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

if __name__ == "__main__":
    condition="osteoporosis"
    Q1_query(condition)

    #condition="rabies"
    #Q1_query(condition)