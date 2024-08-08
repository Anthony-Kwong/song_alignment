#script to make pairwise alignment matrices for all bird pairs, used for pairwise likelihood calculations

library(magrittr)
library(aphid)

source("./functions/align_mat.R")

#read in data
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")

#get all birds
birds = unique(bird_songs$Bird.ID)
#generate all pairs
bp = t(combn(birds, 2, function(x) sort(x), simplify = TRUE))

#loop through pairs and get alignment matrices
for(i in 1:nrow(bp)){
  sam = bp[i,]
  res = align_mat(b1 = sam[1], b2 = sam[2])
  name = paste("./results/bp",i,".RDS", sep="")
  saveRDS(res, file = name)
}