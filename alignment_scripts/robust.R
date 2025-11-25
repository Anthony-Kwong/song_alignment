#Pipeline for comparing different MSA methods on birdsong

#libraries
library("aphid") 
library('magrittr')

#read in data
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
lines = unique(bird_songs$Line)

#load functions
source("./functions/og_order.R")
source("./functions/min_entropy.R")

#pHMMs ---- 

#tuning parameters
lambda = seq(from = 0, to = 100, by = 50)
max_scale = c(1,1.25,1.5) #scale the max number of modules, by the longest length seq
#progress = c(T,F)
tune_params = expand.grid(lambda = lambda, max_scale = max_scale)
#thresh = c(0.4,0.5,0.6)
#tune_params = expand.grid(threshold = thresh, max_scale = max_scale, progress = progress)
#tune_params = expand.grid(threshold = thresh, max_scale = max_scale)

#thresh = seq(from = 0.1, to = 0.9, by = 0.1)

#loop for every song lineage

msa_scores = list()
for(i in 1:length(lines)){
  print(i)
  #filter for birds of the one lineage
  filtered_bird = bird_songs %>%
    dplyr::filter( Line == lines[i])
  
  #get Bird IDs for labeling alignment
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
  
  #get median sequence length (for scaling the max number of modules)
  min_song_len = min(sapply(bird_songs_split,length))
  
  #fit model for different lambda
  
  line_scores = rep(NA, nrow(tune_params))
  for(k in 1:nrow(tune_params)){
    #get set of parameters to fit pHMM
    params = tune_params[k,]
    #fit model
    song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
                            inserts = "map", lambda = params$lambda, maxsize = ceiling(min_song_len*params$max_scale),
                            progressive = FALSE, maxiter = 100)
    
    #threshold method
    # song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
    #                         inserts = "threshold", threshold = tune_params$threshold , maxsize = ceiling(min_song_len*params$max_scale),
    #                         progressive = FALSE, maxiter = 50)
    
    alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
    #retain original order as in bird_songseqs
    A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
    #compute min entropy score
    line_scores[k] = min_entropy(A)
    #save alignment as a fasta file
    alignment_fasta = bio3d::as.fasta(A, id = bIDs)
    fname = paste("./results/fasta/robust_test/pHMM/",lines[i],"lambda_",params$lambda,
                  "_maxscale_",params$max,".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
  }
  #save line results
  msa_scores[[i]] = tibble::tibble(tune_params, score = line_scores)
}

msa_scores = do.call(rbind, msa_scores)
readr::write_csv(msa_scores, file = "./results/msa_scores/phmm_scores.csv")

#20:43 start

#plot the fasta files to check ----
library(ggplot2)
library(ggmsa)
fastas = list.files("./results/fasta/robust_test/pHMM/")
for(i in 1:length(fastas)){
  fname = paste("./results/fasta/robust_test/",fastas[i], sep ="")
  plot = ggmsa(fname, color = "LETTER", seq_name = TRUE, char_width = 0.2) + geom_msaBar() 
  plotname = paste("./results/fasta/robust_test/plot/",fastas[i],".png", sep = "")
  ggsave(plot, file = plotname)
}

#table of song length statistics for every lineage ----
#get song length statistics for every lineage

bird_songs %>%
  dplyr::group_by(Line) %>%
  dplyr::summarise(
    median_length = median(nchar(note.seq), na.rm = TRUE),
    min_length    = min(nchar(note.seq), na.rm = TRUE),
    max_length    = max(nchar(note.seq), na.rm = TRUE),
    .groups = "drop"
  )

#gibbs aligner ----

source("./alignment_scripts/gibbs_align/gibbs_align.R")

#set the number of w values (start from the length of the shortest sequence and then substracting)
w_vals = c(0,1,2)

gibbs_scores = list()
for(i in 1:length(lines)){
  print(i)
  #filter for birds of the one lineage
  filtered_bird = bird_songs %>%
    dplyr::filter( Line == lines[i])
  
  #get Bird IDs for labeling alignment
  IDs = filtered_bird$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = filtered_bird$note.seq
  
  #get minimum song length
  min_len = min(sapply(bird_songsseqs, nchar))
  
  #fit model for different w
  line_scores = rep(NA, length(w_vals))
  for(k in 1:length(w_vals)){
    #fit model
    gibbs.model <- gibbs_align(S = bird_songsseqs, w = min_len - w_vals[k], iter = 100)
    #compute min entropy score
    line_scores[k] = min_entropy(gibbs.model)
    #save alignment as a fasta file
    alignment_fasta = bio3d::as.fasta(gibbs.model, id = bIDs)
    fname = paste("./results/fasta/robust_test/gibbs/",lines[i],"_wminus_",w_vals[k],".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
  }
  #save line results
  gibbs_scores[[i]] = tibble::tibble(wminus = w_vals, score = line_scores, line = lines[i])
}


gibbs_res = do.call(rbind, gibbs_scores)
readr::write_csv(gibbs_res, file = "./results/msa_scores/gibbs_scores.csv")

#progressive alignment ----

source("./alignment_scripts/dynami_program/progressive_align.R")

#set parameters

#tuning parameters
match = c(1,2)
mismatch = c(-1,-2)
gap = c(-1,-2)
#form parameter grid
tune_params = expand.grid(match = match, mismatch = mismatch, gap = gap)

#loop through lineages
dynam_scores = list()
for(i in 1:length(lines)){
  print(i)
  #filter for birds of the one lineage
  filtered_bird = bird_songs %>%
    dplyr::filter( Line == lines[i])
  
  #get Bird IDs for labeling alignment
  IDs = filtered_bird$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = filtered_bird$note.seq
  
  #fit model for different w
  line_scores = rep(NA, nrow(tune_params))
  for(k in 1:nrow(tune_params)){
    params = tune_params[k,]
    #fit model
    dynam_model = progressive_align(S = bird_songsseqs, match = params$match, mismatch = params$mismatch, gap = params$gap)
    #compute min entropy score
    line_scores[k] = min_entropy(dynam_model)
    #save alignment as a fasta file
    alignment_fasta = bio3d::as.fasta(dynam_model)
    fname = paste("./results/fasta/robust_test/dynam_prog/",lines[i],"_match_",params$match
                  ,"_mismatch_", params$mismatch,"_gap_",params$gap ,".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
  }
  #save line results
  dynam_scores[[i]] = tibble::tibble(tune_params, score = line_scores)
}



#old code below----

#investigate robustness of alignments to different parameters

bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
birds = bird_songs$Bird.ID

js = "JS0329"

birds = c("JS0329","JS0037")

js_songs = bird_songs %>%
  dplyr::filter(Bird.ID == js)

setwd("~/Documents/GitHub/song_alignment/results/fasta")

#we are changing match states to insertion states, but this doesnt change the alignment mat A. 

#change maxsize

#affine gap penalties

# 
# 
# #testing parameters----
# 
# #map method for gaps
# lambda = c(0,5,10)
# 
# #go through every bird, generate PHMM alignment with different lambdas and plot
# purrr::walk(birds, function(js){
#   #set model parameters to test----
#   
#   #filter to bird of interest's songs
#   js_songs = bird_songs %>%
#     dplyr::filter(Bird.ID == js)
#   
#   #preprocessing----
#   
#   #get note sequences as long strings
#   bird_songsseqs = js_songs$note.seq
#   
#   #get sequences in split strings
#   bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
#   bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
#   
#   #get the alphabet
#   comb = paste(bird_songsseqs, collapse = "")
#   letters = unique(strsplit(comb, "")[[1]])
#   
#   ###alignment with different parameters
#   
#   # tidyr::crossing(lambda, opt_method)
#   
#   #loop through the 2 optimization procedures
#   
#   #fit the models for different lambdas
#   bw_models = lapply(lambda, function(l){
#     derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
#                inserts = "map", lambda = l)
#   })
#   
#   #create a new directory for the bird if it doesn't exist already
#   bird_dir = paste("./robust_test/map/",js, sep = "")
#   dir.create(bird_dir, showWarnings = FALSE)
#   
#   #loop thru lambdas and generate plots
#   for(i in 1:length(lambda)){
#     alignment = align(bird_songs_split, model = map_models[[i]], seqweights = NULL, residues = letters)
#     A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
#     alignment_fasta = bio3d::as.fasta(A)
#     fname = paste(bird_dir,"/", js,".fasta",sep="")
#     print(fname)
#     bio3d::write.fasta(alignment_fasta, file = fname)
#     
#     #make plot
#     plot = ggmsa(fname, color = "LETTER") + geom_msaBar() 
#     plotname = paste(bird_dir,"/",js,"_",lambda[i],".png", sep = "")
#     ggsave(plot, file = plotname)
#   }
# })
# 
# 
