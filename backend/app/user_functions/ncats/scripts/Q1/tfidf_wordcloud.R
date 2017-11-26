# Word Clouds

# Word cloud of tfidf


library(optparse)
library(wordcloud)

option_list = list(
  make_option(c("-a", "--abstract_1"), type="character", default=NULL, 
              help="Abstract file name", metavar="character"),
  make_option(c("-b", "--abstract_2"), type="character", default=NULL), 
              help="Second abstract file name", metavar="character"),
  make_option(c("-t", "--tfidf"), action="store_true", default=FALSE, 
            help="Abstract file name", metavar="character"),
  # Commonality cloud or Comparison cloud
  # tfidf
  
  make_option()
make_option(c("-p", "--prefix"), type="character", default="abstracts", 
            help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

file_a = opt$abstract_1
file_b = opt$abstract_2

tfidf = opt$tfidf

prefix = opt$prefix

file_a = '/Users/gmcinnes/src/ncats/data/data_dir/tfidf.hiv_infections.txt'
file_b = '/Users/gmcinnes/src/ncats/data/data_dir/abstracts.sickle_cell_trait.txt'

a_terms <- read.table(file_a, header=T, stringsAsFactors=F, quote=NULL, row.names=NULL)
#a_abstracts <- read.table(file_a, header=T, stringsAsFactors=F, quote=NULL, row.names=NULL)
##a_abstracts <- a_abstracts[['ABSTRACT']]
#a_terms <- format_abstracts(a_abstracts, "query_a")

b_terms <- NULL
if (!is.null(file_b)) {
  #b_abstracts <- read.table(file_b, header=T, stringsAsFactors=F, quote=NULL)
  #b_abstracts <- b_abstracts[['ABSTRACT']]
  #b_terms <- format_abstracts(b_abstracts, "query_b")
  b_terms <- read.table(file_b, header=T, stringsAsFactors=F, quote=NULL, row.names=NULL)
}


# Just one
if (is.true(tfidf)) {
  wordcloud(word=a_terms$word, freq=a_terms$tf_idf, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
} else {
  wordcloud(word=a_terms$word, freq=a_terms$n, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
}







# TF IDF word cloud examples
malaria_words = disease_terms[disease_terms$disease == "malaria",] 
wordcloud(word=malaria_words$word, freq=malaria_words$n, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
wordcloud(word=malaria_words$word, freq=malaria_words$tf_idf, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))



sc_words = disease_terms[disease_terms$disease == "sickle_cell_disease",] 
wordcloud(word=sc_words$word, freq=sc_words$n, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
wordcloud(word=sc_words$word, freq=sc_words$tf_idf, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))



sc_words = disease_terms[disease_terms$disease == "sickle_cell_disease",] 
malaria_words = disease_terms[disease_terms$disease == "malaria",] 

head(sc_words)

both = merge(x = sc_words[c("word", "tf_idf")], y = malaria_words[c("word", "tf_idf")], by = "word", all = TRUE)
both[is.na(both)] <- 0
names(both) <- c("word", "sickle_cell", "malaria")
both_m <- as.matrix(both[c("sickle_cell", "malaria")])
rownames(both_m) <- both$word

wordcloud(word=both$word, freq=both$tf_idf, max.words = 100, random.order=FALSE, colors=brewer.pal(8, "Dark2"))



wordcloud(both_m, min.freq =3, scale=c(5, .2), random.order = FALSE, random.color = FALSE, colors= c("indianred1","indianred2","indianred3","indianred"))


comparison.cloud(both_m, max.words=400, random.order=FALSE)

commonality.cloud(both_m, random.order=FALSE, max.words=200, colors=brewer.pal(8, "Dark2"))


