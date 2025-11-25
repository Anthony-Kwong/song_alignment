#benchmark test for alignment (deprecated)

#simulate some sequences
mean_len = 10
sd_len = 2
n = 20
alphabet = 5
motif = "ABBA"
set.seed(1066)
songs = list()
for(i in 1:n){
  #randomise length of new sequence, ensure sequence is at least 4 letters long
  seq_length = max(floor(rnorm(n = 1, mean = mean_len, sd = sd_len)),4)
  base_song = sample(LETTERS[1:alphabet], size = seq_length, replace = T)
  #randomise start position of motif
  pos = sample(2:(seq_length-2), 1)
  #insert motif into base song
  bs1 = paste0(base_song[1:pos], collapse = "")
  #fix bug, string reads backward at edges
  bs2 = paste0(base_song[(pos+1):seq_length], collapse = "")
  fsong = paste0(bs1, motif, bs2)
  songs[[i]] = strsplit(fsong, "")[[1]]
}

#fit phmm
library("aphid")

#fit phmm model

#sort sequences from longest to shortest due to avoid error in align which is called by derivePHMM
ordered_songs <- songs[order(lengths(songs))]
song.PHMM <- derivePHMM(ordered_songs, residues = LETTERS[1:alphabet], pseudocounts = "Laplace")
#plot phmm statespace
plot(song.PHMM)

#simulate sequences for training using fitted phmm
sim <- list(length = 10)
suppressWarnings(RNGversion("3.5.0"))
set.seed(9999)
for(i in 1:10) sim[[i]] <- generate(song.PHMM, size = mean_len)
sim

#we can optimise phmm 
song2.PHMM <- train(song.PHMM, sim, method = "BaumWelch", 
                       deltaLL = 0.01, seqweights = NULL)

#use Viterbi to get path
path = Viterbi(songs2.PHMM, songs[[1]])$path
#rename path in Delete, Match, Insert space
c("D", "M", "I")[path + 1]

align(songs, model = song2.PHMM, seqweights = NULL, residues = LETTERS[1:alphabet])

#compare with gibbs
source("./gibbs_align.R")
x = lapply(songs, function(s){
  paste0(s, collapse = "")
})
gibbs_align(unlist(x), w = 4)

al_res = list(NA)
for(k in 1:100){
  al_res[[k]] = gibbs_align(unlist(x), w = 4)
}

scores = sapply(al_res, function(tab){
  tab[2]
})

best_gibbs = which.max(unlist(scores))
al_res[[best_gibbs]]

#try on song data----
bird_songs = readr::read_csv("./data/NoteSequences.csv")

b37 = bird_songs %>%
  dplyr::filter(bird.num == 37)

b37_songs <- b37$note.seq

#get in form for align function
b37_a = lapply(b37_songs, function(s){strsplit(s, "")[[1]]})
b37_a <- b37_a[order(lengths(b37_a))]


#get the alphabet
comb = paste(b37_songs, collapse = "")
letters = unique(strsplit(comb, "")[[1]])
library(aphid)
song.PHMM <- derivePHMM(b37_songs, residues = letters, pseudocounts = "Laplace")
#get mean song length
mean_len = floor(nchar(gsub("[^A-Za-z]", "", comb))/length(b37_songs))
#plot phmm statespace
plot(song.PHMM)

#simulate sequences for training using fitted phmm
sim <- list(length = 10)
suppressWarnings(RNGversion("3.5.0"))
set.seed(9999)
for(i in 1:10) sim[[i]] <- generate(song.PHMM, size = mean_len)
sim

#we can optimise phmm 
song2.PHMM <- train(song.PHMM, sim, method = "BaumWelch", 
                    deltaLL = 0.01, seqweights = NULL)

plot(song2.PHMM)
alignment_37 = align(b37_a, model = song2.PHMM, seqweights = NULL, residues = letters)


y = as.data.frame(ab37)

library(ggplot2)
alignment_plot <- ggplot(y) +
  geom_text(aes(x = position, y = sequence, label = residue), size = 5) +
  theme_bw() +
  labs(x = "Position", y = "Sequence")


#print example alignment
for(i in 1:nrow(ab37)){
  print(paste0(ab37[i,], collapse = ""))
}


