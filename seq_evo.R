#alignment workflow part 2 (sequence evolution)
#infer transition matrix

library("aphid")
library("ggmsa")
library(magrittr)


bird_songs = readr::read_csv("./data/NoteSequences.csv")
lines = unique(bird_songs$Line)

#fit phmm on a single line ----
#input: name of line as a string
#output: the fitted phmm model

line_align <- function(line){
  filtered_bird = bird_songs %>%
    dplyr::filter( Line == line)
  
  #get note sequences as long strings
  bird_songsseqs = filtered_bird$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  return(song.PHMM)
}

#retrain parameters of a phmm using select sequences only
#input: a PHMM, ID number of select bird to filter songs for
#output: new fitted PHMM

model_refit <- function(model, ID){
  filtered_bird = bird_songs %>%
    dplyr::filter(bird.num == ID)
  
  #get note sequences as long strings
  bird_songsseqs = filtered_bird$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  #retrain PHMM using our select songs
  refit.PHMM <- train(model, y = bird_songs_split, method = "BaumWelch", 
                         deltaLL = 0.01, seqweights = NULL)
  
  return(refit.PHMM)
}

#example using a blue bird 37----
blue_model = line_align(line = "Blue")
b37_model = model_refit(model = blue_model, ID = 37)
#extract emission probabilities
apply(b37_model$E, 2,exp)

#get emission probs for every bird in blue line
x = bird_songs %>%
  dplyr::filter(Line == "Blue")
blue_birds = unique(x$bird.num)

blue_phmms = lapply(blue_birds, function(id){model_refit(ID = id, model = blue_model)})
#extract emission log-probabilities and transform back into regular probability space
blue_pmfs = lapply(blue_phmms, function(x){ apply(x$E, 2, exp) })

pmf1 = rowSums(blue_pmfs[[1]])/ncol(blue_pmfs[[1]])
pmf2 = rowSums(blue_pmfs[[2]])/ncol(blue_pmfs[[2]])
