# Compute semantic similarity between a provided term and genetic diseases downloaded from PubMed
# Altman Lab

# Question 1 genetic disease search
.libPaths()
.libPaths(c(.libPaths(),'/home/admin/anaconda3/envs/iris/lib/R/library/'))
.libPaths(c(.libPaths(), '/home/admin/anaconda3/envs/iris/lib/R/library/', '/home/admin/R/x86_64-pc-linux-gnu-library/3.3', 
	                '/usr/local/lib/R/site-library', '/usr/lib/R/site-library', '/usr/lib/R/library'))
.libPaths()

library(PubMedWordcloud)
library(optparse)
library(lsa)
library(dplyr)
library(tidytext)
library(doParallel)

# Set command line options
option_list = list(
  make_option(c("-q", "--query"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character", default=NULL),
  make_option(c("-p", "--prefix"), type="character", default="temp", 
              help="output file name [default= %default]", metavar="character")
); 

# Parse options
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

query = opt$query

write(paste("Computing semantic similarity for: ", query), stderr())

data_dir = opt$data_dir

# Functions
######################################################
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

# Function to format condition query
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
######################################################

# Get a clean version of the query for saving files
clean_q <- clean_query(query)

# 1. Load the genetic disease term frequencies
genetic_disease_rds <- paste(data_dir, "genetic_disease_tibble.rds", sep="/")
disease_terms <- readRDS(genetic_disease_rds)

# 2. Check if the query is in the genetic disease list
disease_names <- unique(disease_terms$disease)
if (!query %in% disease_names) {
  # 3. If not, get abstracts for the query term?
  # 3a. Check if abstracts are downloaded and in the data directory
  query_abstract_filename <- paste("abstracts", clean_q, "txt", sep=".")
  query_abstract_file_path <- paste(data_dir, query_abstract_filename, sep="/")
  query_abstracts = NA
  
  if (file.exists(query_abstract_file_path)) {
    # Read in the abstracts
    write("Found abstracts locally", stderr())
    query_abstracts <- read.table(query_abstract_file_path, header=T, stringsAsFactors=F, sep="\t", quote=NULL)
    query_abstracts <- query_abstracts[['ABSTRACT']]
  } else {
    # 3b. If not, download them from PubMed  
    write("Downloading abstracts from PubMed", stderr())
    query_abstracts <- getAbstracts(getPMIDsByKeyWords(keys = query, dFrom = 2000, dTo = 2016)[1:1000])
    query_abstracts <- unlist(query_abstracts)
    # 3ba. Save the abstracts to the data directory for future use
    abstracts_df <- data.frame(query_abstracts)
    names(abstracts_df) <- "ABSTRACT"
    abstracts_df <- data.frame(lapply(abstracts_df, function(x) {
      gsub("#", "", x)
    }))
    write.table(abstracts_df, file=query_abstract_file_path, quote=F, row.names = F) 
  }
  # Check if there are no abstracts or if there are too few
  if (length(query_abstracts) < 5) {
    # Should do something to tell python that it broke.  Maybe save a run diagnostics object and load it with ryp2. 
    write("Not enough abstracts found.  Exiting.", stderr())
    quit()
  }
  
  # 3c. Process the abstracts and append them to the genetic disease terms
  query_disease_terms <- format_abstracts(query_abstracts, clean_q)
  disease_terms <- rbind(query_disease_terms, disease_terms)
}


# 4. Calculate tf-idf on the disease terms
disease_terms <- disease_terms %>%
  bind_tf_idf(word, disease, n)
disease_terms

# Save the tfidf values to a files
tfidf_filename <- paste("tfidf", clean_q, "rds", sep=".")
tfidf_file_path <- paste(data_dir, tfidf_filename, sep="/")
saveRDS(disease_terms, file=tfidf_file_path) 


# Extract query terms
query_df <- disease_terms[disease_terms$disease == clean_q,][c("word", "tf_idf")]
names(query_df)[2] <- clean_q

query_df_all_words <- data.frame(unique(disease_terms$word))
names(query_df_all_words) <- "word"
query_df_all_words <- left_join(query_df_all_words, query_df)
query_df_all_words[is.na(query_df_all_words)] <- 0

# Create cluster
cl <- makeCluster(4)
registerDoParallel(cl)
write("Computing cosine similarities", stderr())

# Compute pairwise similarity for each genetic disease against the query disease
similarities <- foreach (i=1:length(disease_names),
              .packages=c('dplyr','lsa'),
              .combine=rbind) %dopar% {
  # Extract disease df
  t <- disease_names[i]
  if (t==clean_q) {
    next()
  }
  t_disease <- disease_terms[disease_terms$disease == t,]
  if(nrow(t_disease) == 0) {
    next()
  }
  t_words <- t_disease[c('word', 'tf_idf')]
  names(t_words)[2] <- t
                
  query_disease_df <- left_join(query_df_all_words, t_words)
  query_disease_df[is.na(query_disease_df)] <- 0
                
  data.frame(disease_a=clean_q,
    disease_b=t, 
    cos_sim=cosine(query_disease_df[[clean_q]], query_disease_df[[t]]), 
    stringsAsFactors=FALSE) 
}  

# Kill the cluster
stopCluster(cl)

# Sort the results
sorted <- similarities %>%
  arrange(desc(cos_sim))

# Save the results to a file
query_sim_file <- paste("similarities", clean_q, "txt", sep=".")
query_sim_file_path <- paste(data_dir, query_sim_file, sep="/")
write.table(sorted, file=query_sim_file_path, quote=F, row.names = F, sep=",") 


