#run dependencies
source("./functions/bck_pdf.R")
source("./functions/string_slice.R")
source("./functions/count_letters.R")
source("./functions/get_alphabet_from_sequences.R")
source("./functions/kmer_bprob.R")
source("./functions/min_entropy.R")
source("./functions/profile.R")
source("./functions/kmer_prob.R")

#' Gibb's aligner algorithm to find a common pattern of width W in a set of input sequences S. 
#'
#' @param S : Vector of strings to align. 
#' @param w : Numeric integer. Width of alignment. 
#' @param iter : Numeric scalar. Number of iterations to carry. 
#'
#' @return : A character matrix of the final alignment. 
#' @export
#'
#' @examples S = c("ABBACC", "CCABBA", "CABBA"), gibbs_align(S, w = 4)
gibbs_align = function(S, w, iter){
  #check inputs ----
  
  #lengths of input sequences
  n = sapply(S, nchar)
  #number of sequences
  k = length(S)
  #ensure that every string in S is at least w long
  if(any(n < w)){
    stop("All input strings must be at least w long.")
  }
  alphabet = get_alphabet_from_sequences(S)
  #compute background probs
  
  #get background probs of A_star
  p = bck_pdf(S)
  
  #initialization ----
  
  #Matrix A contains the current best alignment
  #Start with initial guess alignment for A, which is just the first window of length w for all input strings
  A = matrix(0, nrow = k, ncol = w)
  #fill A with substrings of S (first w elements)
  for(i in 1:k){
    A[i,] = string_slice(S[i], start = 1, end = w)
  }
  #score the alignment 
  A_score = min_entropy(A)
  A_best = A
  
  #gibb's regime 
  
  for(i in 1:iter){
    #predictive update ----
    #print(i)
    #sample 1 row, z to remove
    z = sample(x = k, size = 1)
    #remove z from A (A_star)
    A_star = A[-z,]
    #get pattern probabilities (qi's) in A_star 
    q = profile(A_star, alphabet = alphabet,pseudocount = 1)
    
    
    #sampling step (perturb a_z) ----
    #save row z
    a_z = S[z]
    #vectorise a_z
    a_z = unlist(strsplit(a_z,split = ""))
    #obtain starting indices for valid windows of length w
    z_starts = seq(length(a_z) - w + 1)
    #get candidate kmers
    kmers = lapply(z_starts, function(x){a_z[x:(w+x-1)]})
    
    #compute weights Ax for every kmer
    Ax = rep(NA, length(kmers))
    for(i in 1:length(kmers)){
      #compute prob of kmer under profile of A_star
      qx = kmer_prob(kmer = kmers[[i]], profile = q)
      #compute prob of kmer under random model
      px = kmer_bprob(kmer = kmers[[i]], background = p)
      #compute weight
      Ax[i] = qx/px
    }
    
    #normalise the weights
    Ax = Ax/sum(Ax)
    n_kmers = length(kmers)
    #sample new substring of z based on the weights Ax
    z_sam = sample(x = seq(n_kmers),size = 1 , prob = Ax)
    #get new kmer for sequence z
    z_kmer = unlist(kmers[z_sam])
    #replace zth row of A
    A[z,] = z_kmer
    
    #score the updated guess A
    score_update = min_entropy(A)
    
    if(score_update < A_score){
      #save the updated guess if it scores better than the previous best guess
      A_best = A
    }
  }
  #return best guess alignment
  return(A_best)
}

data = c("ABBACC", "BBABBA","CABBA")
output = gibbs_align(data, w = 4, iter = 500)

data = c("AAABBACC", "CCBBABBA","CABBA")
output = gibbs_align(data, w = 4, iter = 500)

