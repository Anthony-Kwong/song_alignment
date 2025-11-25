#Script for aligning a set of songs using PHMMs. Here we make a separate alignment for every singer and 
#save the fasta files. 

library(aphid)
library("ggmsa")
library('magrittr')

setwd("~/Documents/GitHub/song_alignment/")

source("./functions/og_order.R")

#try on song data----
bird_songs = readr::read_csv("./data/NoteSequences.csv")

bird_songs = bird_songs %>%
  dplyr::mutate(singer = paste("JS",bird.num, sep=""))

birds = unique(bird_songs$singer)

for(i in 1:length(birds)){
  print(i)
  filtered_bird = bird_songs %>%
    dplyr::filter(singer == birds[i])
  
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
  
  #plot PHMM model, for debugging only
  #plot(song.PHMM)
  
  #generate multiple alignment using our fitted PHMM model 
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  
  #retain original order as in bird_songseqs
  A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
  
  #convert into fasta file and save
  alignment_fasta = bio3d::as.fasta(A)
  fname = paste("./data/songs_fasta/",birds[i],".fasta",sep="")
  print(fname)
  bio3d::write.fasta(alignment_fasta, file = fname)
}

library(ggplot2)
library(ggmsa)
#plot a fasta alignment
ggmsa("./data/songs_fasta/JS37.fasta", color = "LETTER") + geom_seqlogo() + geom_msaBar() 

#plot alignments per individual
setwd("./data/songs_fasta/")
fastas = list.files()
for(i in 1:length(fastas)){
  fname = paste("./",fastas[i], sep ="")
  plot = ggmsa(fname, color = "LETTER") + geom_seqlogo() + geom_msaBar() 
  plotname = paste("../alignment_plots/",fname,".png", sep = "")
  ggsave(plot, file = plotname)
}

##main script ends here. 

#following is some code for playing around (to be removed in the final version)



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

#simulate sequences for training using fitted phmm (not really necessary, here for illustrative purposes only)
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

#more examples for playing around

#convert song sequences to fasta (example)
b37_a

alignment_37
a_fasta = bio3d::as.fasta(alignment_37)
bio3d::write.fasta(a_fasta, file = "./data/songs_fasta/JS0037.fasta")
ggmsa("./data/songs_fasta/JS0037.fasta", color = "LETTER") + geom_seqlogo() + geom_msaBar()

library(seqinr)
write.fasta(b37_a, names = paste("JS0037", seq(1:10),sep="_"), file.out = "./data/songs_fasta/JS0037.fasta", as.string = F)

song_sequences <- system.file("./data/songs_fasta/JS0037.fasta")
ggmsa("./data/songs_fasta/JS0037.fasta") + geom_seqlogo() + geom_msaBar()

protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
ggmsa(protein_sequences, start = 221, end = 280, char_width = 0.5, seq_name = T) + geom_seqlogo() + geom_msaBar()
