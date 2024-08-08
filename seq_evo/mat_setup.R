#script for setting out the matrices for seq evo. problem

library(magrittr)
library(aphid)

source("./functions/og_order.R")
source("./functions/col_pseudocounts.R")
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")

#get all father son pairs----
birds = unique(bird_songs$Bird.ID)
sf_pairs = lapply(birds, function(b){
  k = bird_songs %>%
    dplyr::filter(Bird.ID == b)
  tibble::tibble(son = b, dad = unique(k$Social.Father))
})

sf_pairs = do.call(rbind, sf_pairs)
sf_pairs = na.omit(sf_pairs)

#get letters globally
global_letters = LETTERS[1:16]

#add social lineages
meta.data = read.csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Files for Anthony/JavaSparrow_Metadata.csv")

sf_pairs$Line = NA
for(i in 1:nrow(sf_pairs)){
  sam = sf_pairs[i,]
  L1 = meta.data[which(meta.data$Bird.ID == sam$son),]$Line
  L2 = meta.data[which(meta.data$Bird.ID == sam$dad),]$Line
  testthat::expect_equal(L1,L2)
  sf_pairs$Line[i] = L1
}
ß
#compute X,Y matrices ----

pseudo_count = 1
X = list() #dad matrix, probs
Y = list() #son matrix, probs

#matrix of counts
cX = list()
cY = list()

#loop thru the 58 sf pairs
for(i in 1:nrow(sf_pairs)){
  d = sf_pairs[i,] 
  son = d$son
  dad = d$dad
  
  #get set of songs to align ----
  son_songs = bird_songs %>%
    dplyr::filter(Bird.ID == son)
  
  dad_songs = bird_songs %>% 
    dplyr::filter(Bird.ID == dad)
  
  #alignment ----
  #get note sequences as long strings
  bird_songsseqs = c("D" = dad_songs$note.seq, "S" = son_songs$note.seq)
  #indices for fathers and sons
  index = substr(names(bird_songsseqs), 1, 1)
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet for the alignment only
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  letters = letters[order(letters)]
  
  #fit model
  song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
  alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
  #retain original order as in bird_songseqs
  A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
  #save alignment matrix
  fname = paste("./seq_evo/pairwise_fasta/",son,"-",dad, ".csv", sep ="")
  print(fname)
  write.csv(A, file = fname)
  
  #convert to probabilities for matrices X,Y -----
  pre_X = A[which(index == "D"),]
  pre_Y = A[which(index == "S"),]
  #convert columns into probs using pseudocounts 
  X[[i]] = apply(pre_X, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = pseudo_count)
  Y[[i]] = apply(pre_Y, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = pseudo_count)
  
  #get matrix of counts as well
  cX[[i]] = apply(pre_X, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = pseudo_count, pro = FALSE)
  cY[[i]] = apply(pre_Y, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = pseudo_count, pro = FALSE)
}

#check dimensions
a = sapply(X, function(k){ncol(k)})
sum(a)

#matrices with probablities
pro_X = do.call(cbind, X)
pro_Y = do.call(cbind, Y)

#matrices with raw counts

saveRDS(pro_X, file ="./results/seq_evo/sf_mat.rds")
saveRDS(pro_Y, file ="./results/seq_evo/son_mat.rds")


#single case example---- (for debug and code building)
c_X = do.call(cbind, cX)
c_Y = do.call(cbind, cY)

saveRDS(c_X, file ="./results/seq_evo/counts_sf_mat.rds")
saveRDS(c_Y, file ="./results/seq_evo/count_son_mat.rds")

##script complete


son = "JS0040"
dad = "JS0003"

#get set of songs to align ----
son_songs = bird_songs %>%
  dplyr::filter(Bird.ID == son)

dad_songs = bird_songs %>% 
  dplyr::filter(Bird.ID == dad)

#alignment ----
#get note sequences as long strings
bird_songsseqs = c(dad_songs$note.seq, son_songs$note.seq)
#indices for fathers and sons
index = c(rep("D",length(dad_songs$note.seq)), rep("S",length(son_songs$note.seq)))

#get sequences in split strings
bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]

#get the alphabet
comb = paste(bird_songsseqs, collapse = "")
letters = unique(strsplit(comb, "")[[1]])
letters = letters[order(letters)]

#fit model
song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch")
alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
#retain original order as in bird_songseqs
A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
#save alignment as a fasta file
alignment_fasta = bio3d::as.fasta(A)
fname = paste("./seq_evo/pairwise_fasta/pair1.fasta")
print(fname)
bio3d::write.fasta(alignment_fasta, file = fname)

#convert to probabilities for matrices X,Y -----
pre_X = A[which(index == "D"),]
pre_Y = A[which(index == "S"),]
#convert columns into probs using pseudocounts 
X = apply(pre_X, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = 1)
Y = apply(pre_Y, MARGIN = 2, col_pseudocounts, letters = global_letters, pcount = 1)

