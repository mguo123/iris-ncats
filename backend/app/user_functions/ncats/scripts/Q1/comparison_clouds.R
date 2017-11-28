# Generate word clouds to compare two diseases

.libPaths(c(.libPaths(), '/home/admin/anaconda3/envs/iris/lib/R/library/', '/home/admin/R/x86_64-pc-linux-gnu-library/3.3',
            '/usr/local/lib/R/site-library', '/usr/lib/R/site-library', '/usr/lib/R/library'))

library(optparse)
library(wordcloud)

option_list = list(
  make_option(c("-a", "--disease_1"), type="character", default=NULL, 
              help="Disease name", metavar="character"),
  make_option(c("-b", "--disease_2"), type="character", default=NULL, 
            help="Second disease name", metavar="character"),
  make_option(c("-f", "--freq_file"), type="character", default=NULL, 
            help="Frequency file name", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="word_cloud", 
              help="Output prefix", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

clean_query <- function(query) {
  clean_query <- tolower(query)
  clean_query <- gsub(" ", "_", clean_query, fixed=TRUE)
  clean_query <- gsub(",", "", clean_query, fixed=TRUE)
  return(clean_query)
}

# Set the command line options as variables
disease_1 = clean_query(opt$disease_1)
disease_2 = clean_query(opt$disease_2)
prefix = opt$prefix
freq_file <- opt$freq_file

# Read in the frequencies tibble
frequencies <- readRDS(freq_file)

# Extract word frequencies for the desired diseases
a_terms <- frequencies[frequencies$disease == disease_1,]
b_terms <- frequencies[frequencies$disease == disease_2,]

# Merge them together and format as a matrix
both = merge(x = a_terms[c("word", "tf_idf")], y = b_terms[c("word", "tf_idf")], by = "word", all = TRUE)
both[is.na(both)] <- 0
names(both) <- c("word", disease_1, disease_2)
both_m <- as.matrix(both[c(disease_1, disease_2)])
rownames(both_m) <- both$word

# Set the output files
comparison_output_file <- paste(opt$prefix, disease_1, disease_2, "comparison.png", sep=".")
commonality_output_file <- paste(opt$prefix, disease_1, disease_2, "commonality.png", sep=".")

# Create a comparison cloud and save it to a file
png(comparison_output_file, width=12, height=8, units="in", res=300)
comparison.cloud(both_m, max.words=300, random.order=FALSE, colors=c('#cc9837', '#81a558'))
dev.off()

# Create a commonality cloud and save it to a file
png(commonality_output_file, width=12, height=8, units="in", res=300)
layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
par(mar=rep(0, 4))
plot.new()
title = toupper(paste(disease_1, "-", disease_2, "commonality cloud"))
text(x=0.5, y=0.5, title, cex=1.5)
commonality.cloud(both_m, random.order=FALSE, max.words=200, colors='#d45d50')
dev.off()


