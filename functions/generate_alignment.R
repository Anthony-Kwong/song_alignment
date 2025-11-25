#' Generate_alignment function
#' 
#' Generates multiple alignment of birdsongs using PHMMs. 
#'
#' @param songs: Dataframe containing selected rows of NoteSequences.csv contained in the data directory. 
#'
#' @return Fasta file of song alignment
#' @export
#'
#' @examples bird_songs = read.csv(./data/NoteSequences.csv)
#' sub_songs = dplyr::filter(bird_songs, clutch == A)
#' generate_alignment(sub_songs)
generate_alignment <- function(songs){
 
  #get Bird IDs for labeling alignment
  IDs = songs$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = songs$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  #get songs into decreasing order
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  #fit model
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  #retain original order as in bird_songseqs
  A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
  
  #save alignment as a fasta file
  alignment_fasta = bio3d::as.fasta(A, id = bIDs)
  return(alignment_fasta)
#  bio3d::write.fasta(alignment_fasta, file = fname)
}

