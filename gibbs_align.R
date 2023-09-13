#implement gibbs' aligner

#example sequences

data =  c("CCCCCCAAABBAAACCCCC", "DDDDDDAAABBAAADDDDDD", "AAAAAAAAAAAABBAAAAAAAAAAA", "BBBBBBBBAAABBAAABBBBBBB", "CCCCCCAAABBAAACCCCC", "DDDDDDAAABBAAADDDDDD", "AAAAAAAAAAAABBAAAAAAAAAAA", "BBBBBBBBAAABBAAABBBBBBB")
w = 4

#Gibbs aligner----

#run dependencies
source("./functions/bck_pdf.R")
source("./functions/pssm.R")
source("./functions/opmf.R")

#input: 
#t:A vector of character strings, w: integer for the length of motif to discover
#output: Local sequence alignment matrix A, scoring matrix S
gibbs_align = function(t, w){
  #initialization ----
  #index of special sequence
  z = 1
  #lengths of input sequences
  n = sapply(t, nchar)
  #initial special sequence (motif)
  t_star = unlist(strsplit(t[z], split=""))
  n_star = n[z]
  #number of sequences
  k = length(t)
  
  #index of non-special sequences
  index = rep(NA,k)
  
  #matrix to hold best current estimate
  A = matrix(0, nrow = k-1, ncol = w)
  
  for(j in 2:k){
    index[j-1] = j
    #guess starting offset
    o = sample(seq(n[j]-w), 1)
    #update with new candidate alignment after applying offset, we need to split this into its separate letters
    new_motif = substr(t[j], start = o+1 , stop = o+w)
    letters = strsplit(new_motif, split = "")
    letter_vec = unlist(letters)
    #apply update
    A[j-1, ] = letter_vec 
  }
  
  #compute background pdf of letters
  bpdf = bck_pdf(t)
  #compute propensity matrix for A, P
  P = pssm(A, bpdf)
  
  #initialise P for next step
  P2 = matrix(0, nrow = nrow(P), ncol = w)
  deltaP = 900
  
  #keep searching until P converges
  while( deltaP > .Machine$double.xmin){
    #initialise offset pmf
    pmf = tibble::tibble(o = c(0, seq(n_star-w)), prob = NA)
    #get pmf of offset values o
    density = sapply(c(0,seq(n_star-w)), function(o){opmf(o, P, t_star, shift = n_star-w)})
    pmf$prob = density
    #select new offset
    ostar = sample(pmf$o, size = 1, prob = pmf$prob )
    #select new special sequence (motif)
    r = sample(seq(k-1), size = 1)
    #add sequence to A
    A[r,] = t_star[(ostar+1):(ostar+w)]
    #store ptr to t_star in index
    y = index[r]
    index[r] = z
    z = y
    #initialise new t_star
    t_star = unlist(strsplit(t[z], split=""))
    n_star = length(t_star)
    P2 = pssm(A, bpdf)
    
    #add perturbance of P for next iteration
    deltaP = abs(sum(P2-P))
    print(deltaP)
    #update P for next iteration
    P = P2
  }
  
  #compute final A
  final_A = rbind(t_star[(ostar+1):(ostar+w)], A)
  
  #compute score
  rows = nrow(final_A)
  #compute scores for each row in final_A
  scores = sapply(seq(rows), function(i){
    s = rep(NA,w)
    al = final_A[i,]
    for(i in 1:w){
      l = al[i]
      p_index = which(rownames(P)==l)
      s[i] = log(P[p_index,i]) #check the log trick that mark suggested later
    }
    rowscore = sum(s)
    return(rowscore)
  })
  #total score
  final_score = sum(scores)
  
  res = list(A= final_A, score = final_score)
  
  return(res)
}

# gibbs_align(data, w = 8)
# 
# al_res = list(NA)
# for(k in 1:100){
#   al_res[[k]] = gibbs_align(data, w = 8)
# }
# 
# scores = sapply(al_res, function(tab){
#   tab[2]
# })


#if the motif is too short and there's not that many sequences, it struggles
