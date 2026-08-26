#Title : comparison pipeline for alignment techniques ----

#libraries
library("aphid") 
library('magrittr')

#read in data
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
lines = unique(bird_songs$Line)

#load functions
source("./functions/og_order.R")
source("./functions/min_entropy.R")

#home-made alignment functions
source("./alignment_scripts/dynami_program/progressive_align.R")
source("./alignment_scripts/gibbs_align/gibbs_aligner_global.R")


#intra-lineage results. Align entire lineages seperately -----

##pHMMs ---- 

#tuning parameters
lambda = seq(from = 0, to = 100, by = 50)
max_scale = c(1,1.25,1.5) #scale the max number of modules, by the longest length seq
tune_params = expand.grid(lambda = lambda, max_scale = max_scale)

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
  
  #get minimum sequence length (for scaling the max number of modules)
  min_song_len = min(sapply(bird_songs_split,length))
  
  score_stats = list()
  for(k in 1:nrow(tune_params)){
    #get set of parameters to fit pHMM
    params = tune_params[k,]
    #fit model
    song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
                            inserts = "map", lambda = params$lambda, maxsize = ceiling(min_song_len*params$max_scale),
                            progressive = FALSE, maxiter = 100)
    
    alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
    #retain original order as in bird_songseqs
    A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
    #compute min entropy score
    line_scores = min_entropy(A)
    #get length
    align_len = ncol(A)
    #put stats together 
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
    #save alignment as a fasta file
    alignment_fasta = bio3d::as.fasta(A, id = bIDs)
    fname = paste("./results/fasta/robust_test/pHMM/intra/",lines[i],"lambda_",params$lambda,
                  "_maxscale_",params$max,".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
  }
  #save line results
  res = tibble::tibble(tune_params, comp = lines[i])
  scores = do.call(rbind, score_stats)
  msa_scores[[i]] = cbind(res, scores)
}

msa_scores = do.call(rbind, msa_scores)
readr::write_csv(msa_scores, file = "./results/msa_scores/phmm_scores.csv")

##gibbs ----

#tuning params
w_minus = c(0,1,2,3,4,5)

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
  score_stats = list()
  for(k in 1:length(w_minus)){
    #fit model
    gibbs.model <- gibbs_aligner_global(S = bird_songsseqs, w = min_len - w_minus[k], iter = 100)
    
    #gather scores data
    line_scores = min_entropy(gibbs.model)
    align_len = ncol(gibbs.model)
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
    
    #save alignment as a fasta file for visualisation
    # alignment_fasta = bio3d::as.fasta(gibbs.model, id = bIDs)
    # fname = paste("./results/fasta/robust_test/gibbs/intra/",lines[i],"_wminus_",w_minus[k],".fasta",sep="")
    # print(fname)
    # bio3d::write.fasta(alignment_fasta, file = fname)
  }
  #gather score_stats
  scores = do.call(rbind, score_stats)
  res = tibble::tibble(wminus = w_minus, line = lines[i])
  
  #save line results
  gibbs_scores[[i]] = cbind(res, scores)
}
gibbs_res = do.call(rbind, gibbs_scores)
readr::write_csv(gibbs_res, file = "./results/msa_scores/gibbs_scores.csv")

##progressive align----

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
  score_stats = list()
  for(k in 1:nrow(tune_params)){
    params = tune_params[k,]
    #fit model
    dynam_model = progressive_align(S = bird_songsseqs, match = params$match, mismatch = params$mismatch, gap = params$gap)
    
    #reordering alignment mat for visualisation
    
    #put strings back into original order
    A = og_order(align_mat = dynam_model, song_seqs = bird_songsseqs)
    #add back sequence names
    rownames(A) = bIDs
    
    #compute min entropy score
    line_scores = min_entropy(A)
    align_len = ncol(A)
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
    
    
    #save alignment as a fasta file
    alignment_fasta = bio3d::as.fasta(A)
    fname = paste("./results/fasta/robust_test/dynam_prog/intra/",lines[i],"_match_",params$match
                  ,"_mismatch_", params$mismatch,"_gap_",params$gap ,".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
  }
  
  #gather score_stats
  scores = do.call(rbind, score_stats)
  res = tibble::tibble(tune_params, line = lines[i])
  
  #save line results
  dynam_scores[[i]] = cbind(res, scores)
}


final_dp_res = do.call(rbind, dynam_scores)
readr::write_csv(final_dp_res, file = "./results/msa_scores/dp_scores.csv")

#inter-lineage scores ----

##Helper functions ----

source("./alignment_scripts/sample_songs.R")

#one round of lineage pair sampling
paired_data = all_lineage_pairs(bird_songs, line_col = "Line")

#number of sampling iterations per lineage pair
set.seed(1066)
nsim = 5
paired_data <- unlist(
  lapply(1:nsim, function(i) {
    all_lineage_pairs(bird_songs, line_col = "Line")
  }),
  recursive = FALSE
)

#set save to T for saving fastas files
save_fasta = F

##phmm ----

#tuning parameters
lambda = seq(from = 0, to = 100, by = 50)
max_scale = c(1,1.25,1.5) #scale the max number of modules, by the longest length seq
phmm_tune_params = expand.grid(lambda = lambda, max_scale = max_scale)

#loop for every pair and generate an alignment
inter_phmm_scores = list()
for(i in 1:length(paired_data)){
  print(i)
  #get data for the run
  pair_dat = paired_data[[i]]
  #get Bird IDs for labeling alignment
  IDs = pair_dat$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  #get note sequences as long strings
  bird_songsseqs = pair_dat$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  #get minimum sequence length (for scaling the max number of modules)
  min_song_len = min(sapply(bird_songs_split,length))
  
  #fit model for different lambda
  score_stats = list()
  for(k in 1:nrow(phmm_tune_params)){
    #get set of parameters to fit pHMM
    params = phmm_tune_params[k,]
    #fit model
    song.PHMM <- derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
                            inserts = "map", lambda = params$lambda, maxsize = ceiling(min_song_len*params$max_scale),
                            progressive = FALSE, maxiter = 100)
    
    alignment = align(bird_songs_split, model = song.PHMM, seqweights = NULL, residues = letters)
    #retain original order as in bird_songseqs
    A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)

    #compute min entropy score
    line_scores = min_entropy(A)
    #get length
    align_len = ncol(A)
    #put stats together 
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
    
    # #option to save fastas
    # if(save_fasta == T){
    #   #save alignment as a fasta file
    #   alignment_fasta = bio3d::as.fasta(A, id = bIDs)
    #   pair_name = names(paired_data)[i]
    #   fname = paste("./results/fasta/robust_test/pHMM/inter/",pair_name,"lambda_",params$lambda,
    #                 "_maxscale_",params$max,".fasta",sep="")
    #   print(fname)
    #   bio3d::write.fasta(alignment_fasta, file = fname)
    # }
  }
  #collect score stats
  scores = do.call(rbind, score_stats)
  pair_name = names(paired_data)[i]
  res = tibble::tibble(phmm_tune_params, pair = pair_name)
  #save results over all pairs
  inter_phmm_scores[[i]] = cbind(res, scores)
}
#save output
inter_phmm_scores = do.call(rbind, inter_phmm_scores)
readr::write_csv(inter_phmm_scores, file = "./results/msa_scores/inter_phmm_scores.csv")

##gibbs ----

#tuning params
w_minus = c(0,1,2,3,4,5)

gibbs_scores = list()
for(i in 1:length(paired_data)){
  print(i)
  
  pair_dat = paired_data[[i]]
  
  #get Bird IDs for labeling alignment
  IDs = pair_dat$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = pair_dat$note.seq
  
  #get minimum song length
  min_len = min(sapply(bird_songsseqs, nchar))
  
  #fit model for different w
  score_stats = list()
  for(k in 1:length(w_minus)){
    #fit model
    gibbs.model <- gibbs_aligner_global(S = bird_songsseqs, w = min_len - w_minus[k], iter = 100)
    
    #gather scores data
    line_scores = min_entropy(gibbs.model)
    align_len = ncol(gibbs.model)
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
  }
  #gather score_stats
  scores = do.call(rbind, score_stats)
  pair_name = names(paired_data)[i]
  res = tibble::tibble(wminus = w_minus,pair = pair_name)
  
  #save line results
  gibbs_scores[[i]] = cbind(res, scores)
}
gibbs_res = do.call(rbind, gibbs_scores)
readr::write_csv(gibbs_res, file = "./results/msa_scores/inter_gibbs_scores.csv")

##progressive align----

#set parameters

#tuning parameters
match = c(1,2)
mismatch = c(-1,-2)
gap = c(-1,-2)
#form parameter grid
pa_tune_params = expand.grid(match = match, mismatch = mismatch, gap = gap)

#loop through lineage pair samples
dynam_scores = list()
for(i in 1:length(paired_data)){
  print(i)
  
  pair_dat = paired_data[[i]]

  #get Bird IDs for labeling alignment
  IDs = pair_dat$Bird.ID
  bIDs = paste0(IDs,"_", ave(IDs, IDs, FUN = seq_along))
  
  #get note sequences as long strings
  bird_songsseqs = pair_dat$note.seq
  
  #fit model for different w
  score_stats = list()
  for(k in 1:nrow(pa_tune_params)){
    params = pa_tune_params[k,]
    #fit model
    dynam_model = progressive_align(S = bird_songsseqs, match = params$match, mismatch = params$mismatch, gap = params$gap)
    
    #reordering alignment mat for visualisation
    
    #put strings back into original order
    A = og_order(align_mat = dynam_model, song_seqs = bird_songsseqs)
    #add back sequence names
    rownames(A) = bIDs
    
    #compute min entropy score
    line_scores = min_entropy(A)
    align_len = ncol(A)
    score_stats[[k]] = tibble::tibble(score = line_scores, length = align_len, norm_score = line_scores/align_len)
  }
  
  #gather score_stats
  scores = do.call(rbind, score_stats)
  pair_name = names(paired_data)[i]
  res = tibble::tibble(pa_tune_params,pair = pair_name)
  #save line results
  dynam_scores[[i]] = cbind(res, scores)
}

final_dp_res = do.call(rbind, dynam_scores)
readr::write_csv(final_dp_res, file = "./results/msa_scores/inter_dp_scores.csv")

