#example 

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
  
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  #JS220 (69) only sings one note type so we can't align it
  
  #plot(song.PHMM)
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  
  alignment_fasta = bio3d::as.fasta(alignment)
  fname = paste("./data/songs_fasta/",birds[i],".fasta",sep="")
  print(fname)
  bio3d::write.fasta(alignment_fasta, file = fname)
}

library(ggplot2)
ggmsa("./data/songs_fasta/JS37.fasta", color = "LETTER") + geom_seqlogo() + geom_msaBar() 

setwd("./data/songs_fasta/")
fastas = list.files()
for(i in 1:length(fastas)){
  fname = paste("./",fastas[i], sep ="")
  plot = ggmsa(fname, color = "LETTER") + geom_seqlogo() + geom_msaBar() 
  plotname = paste("../alignment_plots/",fname,".png", sep = "")
  ggsave(plot, file = plotname)
}



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
