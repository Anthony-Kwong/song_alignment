#' align_mat
#' 
#' Aligns the recordings of a pair of birds. Outputs matrix where the columns were aligned positions 
#' and the rows are the note classes (A - P). Entries are raw counts of the note classes for every 
#' position.
#'
#' @param b1 : BirdID as string. 
#' @param b2 : BirdID as string.
#'
#' @return : List of 2 matrices, one for every bird. Columns are the alignment positions. Rows are the 16 
#' note classes. Entries are the counts. 
#' @export
#'
#' @examples align_mat("JS0003","JS0087")

#load in bird data
library(magrittr)
library(aphid)
source("./functions/og_order.R")
source("./functions/col_pseudocounts.R")
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")

#main function starts here

align_mat <- function(b1,b2, doTests = F){
  
  #set letter range, this was based on the whole population which had 16 note classes
  global_letters = LETTERS[1:16]
  
  #get set of songs to align ----
  b1_songs = bird_songs %>%
    dplyr::filter(Bird.ID == b1)
  
  b2_songs = bird_songs %>% 
    dplyr::filter(Bird.ID == b2)
  
  #alignment ----
  #get note sequences as long strings
  bird_songsseqs = c("X" = b1_songs$note.seq, "Y" = b2_songs$note.seq)
  #indices for fathers and sons
  index = substr(names(bird_songsseqs), 1, 1)
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  #sort in descending order because that is what derivePHMM likes
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet for the alignment only, to use with derivePHMM
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  letters = letters[order(letters)]
  
  #fit model
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  #retain original order as in bird_songseqs
  A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
  
  pre_b1 = A[which(index == "X"),]
  pre_b2 = A[which(index == "Y"),]
  
  b1.mat = apply(pre_b1, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = 0, pro = FALSE)
  b2.mat = apply(pre_b2, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = 0, pro = FALSE)
  
  #alternate output if we are doing tests ----
  if(doTests){
    test_obj = list(bird_songs_split, letters)
    return(test_obj)
  }
  
  #continue with the function ----
  
  #store the alignment matrices
  alignment = list(pre_b1, pre_b2)
  #store the note counts corresponding to the alignment matrices
  counts = list(b1.mat, b2.mat)
  #add bird IDs
  names(alignment) = c(b1,b2)
  names(counts) = c(b1,b2)
  #put everything together
  output = list(alignment = alignment, counts = counts)
  return(output)
}

#tests----

b1 = "JS0003"
b2 = "JS0087"

output = align_mat(b1,b2, doTests = TRUE)
#check correct songs are being passed ----

#filter for the bird
b_data = lapply(c(b1,b2), function(b){
  bird_songs %>% dplyr::filter(Bird.ID == b)
})
#get the note sequences
b_songs = lapply(b_data, function(d){
  data = d
  data$note.seq
})
song_seqs = c("X" = b_songs[[1]], "Y" = b_songs[[2]])
#split strings into characters
song_seqs_sp = lapply(song_seqs, function(s){strsplit(s, "")[[1]]})
#order songs in descending order in length
ordered_songs = song_seqs_sp[order(lengths(song_seqs_sp), decreasing = T)]
testthat::expect_equal(ordered_songs, output[[1]])

#check correct letters are used during model fitting
s = unlist(b_songs, recursive = F)
p = paste(s, collapse = "")
output_letters = unique(strsplit(p, "")[[1]])
L = output_letters[order(output_letters)]
testthat::expect_equal(L, output[[2]])
