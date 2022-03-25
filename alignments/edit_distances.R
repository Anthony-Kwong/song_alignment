pacman::p_load(stringdist,magrittr)
#edit distance experiments

#run script to get sequences
source("./alignment.R")

#dataframe of song sequences
seq.df

#get small sample of songs for starters
ex_songs = seq.df$note.seq[1:10]

#pairwise edit distances
dist_mat <- stringdistmatrix(ex_songs, ex_songs, 
                             #use Levenshtein distance, minimum number of single character edits to transform strings
                             method="lv")

#we can look at within individual distance, take the mean of intra distances
birds = unique(seq.df$bird.num)

intra_dist = sapply(birds, function(ID){
  #filter data by bird
  fil_data = seq.df %>%
    dplyr::filter(bird.num == ID)
  #compute distances between the songs from each individual
  dist_mat <- stringdistmatrix(fil_data$note.seq, fil_data$note.seq, method = "lv")
  #calculate mean distances
  mean_dist <- mean(dist_mat[lower.tri(dist_mat)])
})

intra_df <- tibble::tibble(bird_ID = as.character(birds), intra_dist)

ggplot2::ggplot(data = intra_df, aes(y = intra_dist, x = bird_ID)) +
  geom_bar(stat = "identity")



#distances between individuals



library(alakazam)

pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
             #need an input distance matrix to show closeness of different individual elements
             dist_mat=getDNAMatrix(gap=0))

#Try out some classic alignment algorithms (e.g. Needleman Wunsch)

#check by hand

head(seq.df)
dist(seq.df)

#Now we are ready to do some alignments

#examples seqs
ex_songs = seq.df$note.seq[1:5]
pairwiseDist(ex_songs)

##edit distances
library(stringdist)
stringdist(ex_songs[1], ex_songs[2], method = 'lv')


d <- stringdistmatrix(c('foo','bar','boo','baz'))
# try 
plot(hclust(d))

stringdistmatrix(c("foo","bar","boo"),c("foo","bar","boo"))