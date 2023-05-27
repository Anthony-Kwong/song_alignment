pacman::p_load(stringdist,magrittr,ggplot2)
#edit distance experiments

#dataframe of song sequences
seq.df = read.csv("./data/NoteSequences.csv")

#get small sample of songs for starters
ex_songs = seq.df$note.seq[1:10]

#exploratory analysis----

#pairwise edit distances between sequences
dist_mat <- stringdistmatrix(ex_songs, ex_songs, 
                             #use Levenshtein distance, minimum number of single character edits to transform strings
                             method="lv")

#we can look at within individual distance, take the mean of intra distances
birds = unique(seq.df$bird.num)

#compute mean of the pairwise differences within recordings of the same bird
intra_dist = sapply(birds, function(ID){
  #filter data by bird
  fil_data = seq.df %>%
    dplyr::filter(bird.num == ID)
  #compute distances between the songs from each individual, distance is levenstein
  dist_mat <- stringdistmatrix(fil_data$note.seq, fil_data$note.seq, method = "lv")
  #calculate mean distances
  mean_dist <- mean(dist_mat[lower.tri(dist_mat)])
})

intra_df <- tibble::tibble(bird_ID = as.character(birds), intra_dist)

#plot of mean Levenstein distances within each individual
#note some birds are much more variable than others
ggplot2::ggplot(data = intra_df, aes(y = intra_dist, x = bird_ID)) +
  geom_bar(stat = "identity")


#edit distances using adist----
match_info = adist(ex_songs, counts = TRUE)

#trafos gives operational steps for get from one string to another
match_info[2]
ops <- attributes(match_info)$trafos

#create a function to compute likelihood using the number of steps D,S,I

#having a play rn, will make this a fct on a separate script ----

op = "SMMSMMMIIIIMMMIMIMMMIMMMMMMMMMMM"
res = sapply(c("D","S","I"), function(x){stringr::str_count(op,x)})
Nd = res[1]
Ns = res[2]
Ni = res[3]
lambda = c(1,1,1) #placeholder for the lambda rates for each operation
L = stringr::str_count(op) #placeholder length for now. Need to sort out which of the 2 lengths to use
pair_loglik = Nd*log(lambda[1]) + Ns*log(lambda[2]) + Ni*log(lambda[3]) + 
  (Nd+Ns+Ni)*log(L) - 
  #issue with logging 0, what if there were no D's?
  log(Nd) - log(Ns) - log(Ni) - 
  (sum(lambda))*L

#need to check model theory with mark, step by step33


#distances between individuals (alakazam method)----

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