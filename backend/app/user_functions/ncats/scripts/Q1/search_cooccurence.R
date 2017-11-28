# Provided with two disease names, 
# 1. search pubmed for coccurences
# 2. Compute tf-idf using the individual diseases as background
# 3. Generate word clouds for the tf-idf weighted terms


.libPaths(c(.libPaths(), '/home/admin/anaconda3/envs/iris/lib/R/library/', '/home/admin/R/x86_64-pc-linux-gnu-library/3.3', 
            '/usr/local/lib/R/site-library', '/usr/lib/R/site-library', '/usr/lib/R/library'))


library(PubMedWordcloud)
library(optparse)
library(lsa)
library(dplyr)
library(tidytext)
library(wordcloud)
library(tm)

option_list = list(
  make_option(c("-a", "--disease_1"), type="character", default=NULL, 
              help="Disease name", metavar="character"),
  make_option(c("-b", "--disease_2"), type="character", default=NULL, 
              help="Second disease name", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character", default=NULL),
  make_option(c("-p", "--prefix"), type="character", default="cooccurence", 
              help="Output prefix", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_dir = opt$data_dir


disease_1 = opt$disease_1
disease_2 = opt$disease_2


# Functions
###################################
# Function to format the query for file names

clean_query <- function(query) {
  clean_query <- tolower(query)
  clean_query <- gsub(" ", "_", clean_query, fixed=TRUE)
  clean_query <- gsub(",", "", clean_query, fixed=TRUE)
  return(clean_query)
}


# Convert abstracts into tidy objects
tidy_up_text <- function(a) {
  text_df <- data_frame(line=1:length(a), text = a)
  tidy_terms<- text_df %>%
    unnest_tokens(word, text)
  data(stop_words)
  tidy_terms <- tidy_terms %>%
    anti_join(stop_words)
  return(tidy_terms)
}

# Function to format abstracts downloaded from pubmed
format_abstracts <- function(abstracts, disease) {
  query_terms <- tidy_up_text(abstracts)
  query_terms$disease = disease
  
  # Get counts for each word by disease
  query_disease_terms <- query_terms %>%
    count(disease, word, sort = TRUE) %>%
    ungroup()
  
  # Get sums for each disease
  total_words <- query_disease_terms %>% 
    group_by(disease) %>% 
    summarize(total = sum(n))
  
  # Add the total count to the disease words
  query_disease_terms <- left_join(query_disease_terms, total_words)
  
  return(query_disease_terms)
}


# Fetch abstracts from PubMed
fetch_abstracts <- function(disease_name) {
  query_abstracts = NA
  clean_q <- clean_query(disease_name)
  query_abstract_filename <- paste("abstracts", clean_q, "txt", sep=".")
  query_abstract_file_path <- paste(data_dir, query_abstract_filename, sep="/")
  if (file.exists(query_abstract_file_path)) {
    # Read in the abstracts
    print("Found abstracts locally")
    query_abstracts <- read.table(query_abstract_file_path, header=T, stringsAsFactors=F, sep="\t", quote=NULL)
    query_abstracts <- query_abstracts[['ABSTRACT']]
  } else {
    # 3b. If not, download them from PubChem  
    print("Downloading abstracts")
    query_abstracts <- getAbstracts(getPMIDsByKeyWords(keys = disease_name, dFrom = 2000, dTo = 2016)[1:1000])
    query_abstracts <- unlist(query_abstracts)

    # Save the abstracts to the data directory
    abstracts_df <- data.frame(query_abstracts)
    print(head(abstracts_df))
    names(abstracts_df) <- "ABSTRACT"
    abstracts_df <- data.frame(lapply(abstracts_df, function(x) {
      gsub("#", "", x)
    }))
    write.table(abstracts_df, file=query_abstract_file_path, quote=F, row.names = F) 
  }
  return(query_abstracts)
}


# Get a term dataframe for a given disease
fetch_terms <- function(disease_name) {
  clean_q <- clean_query(disease_name)
  query_abstracts = fetch_abstracts(disease_name)
  # Check if there are no abstracts or if there are too few
  if (length(query_abstracts) < 10) {
    # Should do something to tell python that it broke.  Maybe save a run diagnostics object and load it with ryp2. 
    print(query_abstracts)
    quit()
  }

  # Process the abstracts and append them to the genetic disease terms
  query_disease_terms <- format_abstracts(query_abstracts, clean_q)
  return(query_disease_terms)
}
###################################################

# Get disease terms for each disease
disease_1_terms <- fetch_terms(disease_1)
disease_2_terms <- fetch_terms(disease_2)

# Get disease terms for cooccurence of terms
both_terms <- fetch_terms(paste(disease_1, disease_2))

# Bind them all together into a single tibble
disease_terms <- rbind(disease_1_terms, disease_2_terms)
disease_terms <- rbind(disease_terms, both_terms)

# Calculate tf-idf on the disease terms
disease_terms <- disease_terms %>%
  bind_tf_idf(word, disease, n)

# Print some out just to see something on the command line
disease_terms[with(disease_terms, order(-tf_idf)), ]

# Get a data frame for the cooccurence of the terms
both_clean <- clean_query(paste(disease_1, disease_2))
both <- disease_terms[disease_terms$disease == both_clean,]

# Set the output files
tfidf_output_file <- paste(opt$prefix, clean_query(disease_1), clean_query(disease_2), "cooccurence_tfidf.png", sep=".")
freq_output_file <- paste(opt$prefix, clean_query(disease_1), clean_query(disease_2), "cooccurence_frequency.png", sep=".")

# Generate and save the word clouds
png(tfidf_output_file, width=12, height=8, units="in", res=300)
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mar=rep(0, 4))
plot.new()
title = toupper(paste(disease_1, "-", disease_2, "abstract cooccurrence cloud - tf-idf"))
text(x=0.5, y=0.5, title, cex=1.5)
wordcloud(word=both$word, freq=both$tf_idf, max.words = 200, random.order=FALSE, colors='#ae8a33')
dev.off()


png(freq_output_file, width=12, height=8, units="in", res=300)
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mar=rep(0, 4))
plot.new()
title = toupper(paste(disease_1, "-", disease_2, "abstract cooccurrence cloud"))
text(x=0.5, y=0.5, title, cex=1.5)
wordcloud(word=both$word, freq=both$n, max.words = 200, random.order=FALSE, colors='#229c81')
dev.off()


# Part 2 - Search for sentences with both terms
# Fetch abstracts containing both diseases
both_abstracts <- fetch_abstracts(paste(disease_1, disease_2))
text_df <- data_frame(line=1:length(both_abstracts), text = both_abstracts)

# Create a term tibble with sentences
tidy_abstract_sentences <- text_df %>%
  unnest_tokens(word, text, token="sentences")

# Subset the tible for sentences containing both diseases
d1 <- tidy_abstract_sentences[grepl(disease_1,tidy_abstract_sentences$word, perl=T, ignore.case = T), ]
d2 <- d1[grepl(disease_2,d1$word, perl=T, ignore.case = T), ]

# Write the sentences to a file
sentence_output_file <- paste(opt$prefix, clean_query(disease_1), clean_query(disease_2), "cooccurence_sentences.txt", sep=".")
write.table(d2, file=sentence_output_file, sep="\t", row.names = F)
