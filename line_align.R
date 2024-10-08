#generate multiple alignments for every song lineage

library("aphid")
library("ggmsa")
library('magrittr')

setwd("~/Documents/GitHub/song_alignment/")

source("./functions/og_order.R")


bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
lines = unique(bird_songs$Line)

#check number of recordings for every bird, for data purposes
bird = unique(bird_songs$Bird.ID)
k = sapply(bird, function(b){
  d = bird_songs %>%
    dplyr::filter(Bird.ID == b)
  nrow(d)
})

#alignment code starts here----

#loop for every song lineage
for(i in 1:length(lines)){
  print(i)
  #filter for birds of the one lineage
  filtered_bird = bird_songs %>%
    dplyr::filter( Line == lines[i])
  
  #get Bird IDs for labelling alignment
  IDs = filtered_bird$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = filtered_bird$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  #fit model
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  #JS220 (69) only sings one note type so we can't align it
  
  #plot PHMM for demonstration only
  #plot(song.PHMM)
  
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  #retain original order as in bird_songseqs
  A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
  
  #save alignment as a fasta file
  alignment_fasta = bio3d::as.fasta(A, id = bIDs)
  fname = paste("./results/fasta/lines_fasta/",lines[i],".fasta",sep="")
  print(fname)
  bio3d::write.fasta(alignment_fasta, file = fname)
}

#setwd("./results/fasta/lines_fasta/")

#plot multiple alignments for every song lineage

library(ggplot2)
fastas = list.files(path = "./results/fasta/lines_fasta/")
for(i in 1:length(fastas)){
  fname = paste("./results/fasta/lines_fasta/",fastas[i], sep ="")
  plot = ggmsa(fname, color = "LETTER", seq_name = TRUE, char_width = 0.2) + geom_msaBar() 
  plotname = paste("./results/alignment_plots/",fastas[i],".png", sep = "")
  ggsave(plot, file = plotname)
}



