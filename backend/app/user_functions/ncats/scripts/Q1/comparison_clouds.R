# Word Clouds

# Given a path to one or two abstract files, generate word clouds


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

disease_1 = clean_query(opt$disease_1)
disease_2 = clean_query(opt$disease_2)

tfidf = opt$tfidf

prefix = opt$prefix

#disease_1 = 'hiv infections'
#disease_2 = 'arts_syndrome'

freq_file <- opt$freq_file
#freq_file <- '/Users/gmcinnes/src/ncats/data/data_dir/tfidf.hiv_infections.rds'

frequencies <- readRDS(freq_file)

a_terms <- frequencies[frequencies$disease == disease_1,]
b_terms <- frequencies[frequencies$disease == disease_2,]

both = merge(x = a_terms[c("word", "tf_idf")], y = b_terms[c("word", "tf_idf")], by = "word", all = TRUE)
both[is.na(both)] <- 0
names(both) <- c("word", disease_1, disease_2)
both_m <- as.matrix(both[c(disease_1, disease_2)])
rownames(both_m) <- both$word

comparison_output_file <- paste(opt$prefix, disease_1, disease_2, "comparison.png", sep=".")
commonality_output_file <- paste(opt$prefix, disease_1, disease_2, "commonality.png", sep=".")

png(comparison_output_file, width=12, height=8, units="in", res=300)
comparison.cloud(both_m, max.words=200, random.order=FALSE)
dev.off()

png(commonality_output_file, width=12, height=8, units="in", res=300)
commonality.cloud(both_m, random.order=FALSE, max.words=200, colors=brewer.pal(8, "Dark2"))
dev.off()



