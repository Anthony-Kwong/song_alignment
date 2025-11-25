#Script for finding motifs using pHMMs (proof of concept)

#Write down manually identified motifs from the Turquoise lineage

#Train pHMM on motifs. Use note repertoire of lineage/bird. 

#Take whole songs and make substrings via sliding window. Put into PHMM and find Viterbi path.

#Compute ratio of Viterbi path score/random model (from background probs).

#Plot ration across song. 

#Do cross validation. 5 fold.

#load libs
pacman::p_load(magrittr, aphid)

#load functions
source("./functions/string_split.R")
source("./functions/string_slice.R")
source("./functions/kmer_bprob.R")
source("./functions/bck_pdf.R")

#record motifs ---- 

#set of manually identified motifs
js119_motifs = c("AAAAAAEBEBEABBE", "AAAAAEBEBEABBE",
                 "AAAAAAEBEBEABBE","AAAAAEBEABBE",
                 "AAAAAAEBEBEABBE", "AAAAAEBEABBE",
                 "AAAAAEBEBEABBE", "AAAAAEBEABBE",
                 "AAAAAAEBEBEABBE", "AAAAAEBEABBE",
                 "AAAAAAEBEBEABBE", "AAAAAE",
                 "AAAAAEBEBEABBE", "AAAAAEBEABBE",
                 "AAAAAAEBEBEABBE", "AAAAAEBEABBE",
                 "AAAAAAEBEABBE", "AAAAAEBEABBE",
                 "AAAAAAEBEBEABBE", "AAAAAAE")

js0300_motifs = c("GBBEIIAA","GBBEIIAA","GBBEIE",
                  "GBBEIIAA","GBBEIIAA", "GBBE",
                  "GBBEIIAA","GBBEIIAA", "GBBEIE",
                  "GBBEIIAA","GBBEIIAA", "GBBE",
                  "GBBEIIAA","GBBEIIAA","GBBEIIAA",
                  "GBBEIIAA","GBBEIIAA",
                  "GBBEIIAA","GBBEIIAA",
                  "GBBEIIAA","GBBEIIAA","GBBEIIAA","GBBE",
                  "GBBEIIAA","GBBEIIAA", "GBBEIIA",
                  "GBBEIIAA", "GBBEIIAA", "GBBEI")

motif_list = list(js119_motifs, js0300_motifs)

#model fitting ----

#read in whole song data set
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")

#birds to test
birds = "JS0300"

#get target note seqs
bird_dat = bird_songs %>%
  dplyr::filter(Bird.ID == birds)

#get songs for testing
songs = bird_dat$note.seq
#combine all strings together
comb = paste(songs, collapse = "")
#get singer's whole alphabet
alphabet = unique(strsplit(comb, "")[[1]])
#get background probabilities
bpdf = bck_pdf(songs)

#process motif strings
motif_train = string_split(js0300_motifs)
#set k parameter (to enable fitting in derivepHMM)
min_len = min(sapply(motif_train, length))
k = 0
if(min_len > 5){
  k = 5
} else {
  k = min_len
}

#train motif pHMM
motif.PHMM <- derivePHMM(motif_train, k= k ,residues = alphabet, pseudocounts = "Laplace", refine = "BaumWelch")

#predict on a set of songs

#loop over all songs
song_scores = list()
for(s in 1:length(songs)){
  #set width of window
  w = motif.PHMM$size
  test_song = unlist(string_split(songs[s]))
  #implement sliding window
  nwins = length(test_song) - w 
  #store scores
  pos_scores = rep(NA, nwins)
  #slide across
  for(i in 1:nwins){
    print(i)
    #slice kmer from test_song
    kmer = string_slice(s = test_song, start = i, end = i + w )
    print(kmer)
    #get score under random model
    random_score = kmer_bprob(kmer = kmer, background = bpdf, log_prob = T)
    #get path
    path = Viterbi(motif.PHMM, kmer)
    pos_scores[i] = path$score-random_score
  }
  song_scores[[s]] = pos_scores
}

#plot song 1 scores
song1 = tibble::tibble(score = song_scores[[1]], pos = seq(length(song_scores[[1]])))

library(ggplot2)
ggplot(song1, aes(x=pos, y=score)) +
  geom_point() +
  geom_vline(xintercept = 9, col = "orangered") +
  geom_vline(xintercept = 19, col = "orangered") +
  geom_vline(xintercept = 27, col = "orangered")

#plot scores for several songs (maybe do CV), plot 5 of the 10. Vertical lines for the staring positions. 
