
library("PubMedWordcloud")
library("RColorBrewer")
library('lsa')


query_terms <- c("Osteoporosis", 
                 "Human Immunodeficiency Virus Infection",
                 "Cholera",
                 "Ebola Virus Infection",
                 "Malaria",
                 "Osteomalacia",
                 "Hypercholesterolemia",
                 "Diabetes Mellitus, Type 2",
                 "Asthma",
                 "Chronic Pancreatitis",
                 "Alzheimer Disease",
                 "Myocardial Infarction",
                 "Duchenne Muscular Dystrophy",
                 "Deficiency of N-glycanase 1",
                 "Alcohol Dependence",
                 "Major Depression",
                 "Niemann Pick Type C",
                 "Huntington Disease",
                 "Alkaptonuria",
                 "Sickle Cell Disease",
                 "Post-Traumatic Stress Disorder",
                 "Cystic fibrosis",  # not included in provided list
)

sorted_terms <- sort(query_terms)

filter_frequent_words <- function(word_df) {
  data(stopwords_en)
  frequent_words <- read.table("/Users/gmcinnes/src/ncats/data/stop_words.txt", header=F, stringsAsFactors = F)
  
  
  return(word_df[!word_df$word %in% c(stopwords_en, frequent_words$V1),])
}

td = tempfile()
dir.create(td)
set.seed(123)
for (i in 1:length(query_terms)) {
  
  term <- query_terms[i]
  print(term)
  pmids = getPMIDsByKeyWords(keys = term, dFrom = 2015, 
                             dTo = 2016)
  # limiting to only 1000 for now
  abstracts = getAbstracts(pmids[1:1000])
  cleanAbs = cleanAbstracts(abstracts)
  # Then again only take the top 1000 words
  #cleanAbs[1:1000,]$word
  cleanAbs[] <- lapply(cleanAbs, as.character)
  cleanAbs <- filter_frequent_words(cleanAbs)
  write( cleanAbs[1:1000,]$word, file=paste(td, term, sep="/"))
}


myMatrix = textmatrix(td, minWordLength=1)
#res <- lsa::cosine(myMatrix[,15], myMatrix[,21])
#res

similarities <- data.frame(disease_a=character(),
                 disease_b=character(), 
                 cos_sim=numeric(), 
                 stringsAsFactors=FALSE) 

for (i in 1:length(sorted_terms)) {
  for (j in 1:length(sorted_terms)) {
    if(i == j) {
      next
    }
    res <- lsa::cosine(myMatrix[,i], myMatrix[,j])
    new_row <- data.frame(disease_a=sorted_terms[i],
                          disease_b=sorted_terms[j], 
                          cos_sim=res, 
                          stringsAsFactors=FALSE) 
    similarities <- rbind(similarities, new_row)
  }
}

ggplot(data = similarities, aes(x = disease_a, y = disease_b)) +
  geom_tile(aes(fill = cos_sim)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors=wes_palette(n=3, name="Zissou"))











pmids = getPMIDsByKeyWords(keys = "sickle cell anemia", dFrom = 2010, 
                           dTo = 2013)
abstracts = getAbstracts(pmids)
cleanAbs = cleanAbstracts(abstracts)
filtered <- filter_frequent_words(cleanAbs)
plotWordCloud(filtered, rot.per = 0, min.freq = 3)


##########################




##########################
setwd('/Users/gmcinnes/src/ncats/bin')
td = tempfile()
dir.create(td)
write( c("dog", "cat", "mouse"), file=paste(td, "D1", sep="/") )
write( c("ham", "mouse", "sushi"), file=paste(td, "D2", sep="/") )
write( c("dog", "pet", "pet"), file=paste(td, "D3", sep="/") )

data(stopwords_en)
myMatrix = textmatrix(td, stopwords=stopwords_en)
myMatrix = lw_logtf(myMatrix) * gw_idf(myMatrix)
myLSAspace = lsa(myMatrix, dims=dimcalc_share())
as.textmatrix(myLSAspace)

##########################
td = tempfile()
dir.create(td)
write( c("HDa","2Pb","2","BxU","BuQ","Bve", "Bve"), file=paste(td, "D1", sep="/"))
write( c("HCK","2Pb","2","09","F","G"), file=paste(td, "D2", sep="/"))

# read files into a document-term matrix



df.team_data <- expand.grid(teams = c("Team A", "Team B", "Team C", "Team D")
                            ,metrics = c("Metric 1", "Metric 2", "Metric 3", "Metric 4", "Metric 5")
)

# add variable: performance
set.seed(41)
df.team_data$performance <- rnorm(nrow(df.team_data))

#inspect
head(df.team_data)





tfidf
pointwise mutual information

Try to check just review articles, or check for more abstracts
mapping rare diseases onto common diseases
we care more abot the rare disease mapping to the genetic disease
