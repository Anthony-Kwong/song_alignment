#' pssm function
#'
#' Compute position specific scoring matrices using an alignment and some background distribution for the letters. Method
#' adapted from Computational Molecular Biology 2021, Durand.
#'
#' @param A: A matrix of letters representing some candidate alignments. The rows are the sequences, columns are the sites.
#' @param bpdf: Table for the background distribution of all the elements. 
#'
#' @return: A numeric propensity matrix representing the likelihood ratios for every site.The rows are the letters, colums are the sites.
#' @export
#'
#' @examples B = matrix(c("A","B"), ncol = 4, nrow = 2)
#' seq = strsplit("ABAB",split="")
#' bpdf = prop.table(table(seq))
#' pssm(B, bpdf)
pssm <- function(A, bpdf){
  #get number of instances of motif
  k = nrow(A)
  #get length of motif strings
  w = ncol(A)
  #get alphabet
  sigma = names(bpdf)
  sigma_size = length(bpdf)
  
  #initialize propensity matrix P
  P = matrix(nrow = sigma_size, ncol = w)
  rownames(P) = sigma
  #initialize pseudocount
  b = 1

  #initialize pseudocounts q----
  q = matrix(nrow = sigma_size, ncol = w)
  
  #initialize denominator of background dists
  p = matrix(nrow = sigma_size, ncol = w)
  
  #loop through columns of A
  for(j in 1:w){
    #get counts of every letter in the column (site)
    site_counts = table(A[,j])
    #loop through all letters
    for(i in 1:sigma_size){
      #print(paste("i ", i, "j ", j))
      #get the letter of interest
      l = rownames(P)[i]
      
      #check if l is in A
      if(l %in% names(site_counts)){
        #get count of l in A
        cij = site_counts[which(names(site_counts)==l)]
      } else {
        #otherwise return count 0
        cij = 0
      }
      
      #add entry to pseudo count matrix
      q[i,j] = (cij + b)/(k -1 + b*sigma_size)
      #add entry to background prob. matrix
      p[i,j] = bpdf[which(names(bpdf)==l)]
    }
  }
  
  #compute propensity matrix P, using q and p
  P = q/p
  rownames(P) = sigma
  return(P)
}

#add test

library(testthat)

#test1
B = matrix(c("A","B"), ncol = 4, nrow = 2)
seq = strsplit("ABAB",split="")
bpdf = prop.table(table(seq))
output = pssm(B, bpdf)

q = matrix((1+1)/(2+1*2), ncol = 4, nrow = 2)
p = matrix(0.5, ncol = 4, nrow = 2)
ans = q/p
rownames(ans) = c("A","B")
expect_equal(output, ans)


#test2
C = matrix(ncol = 5, nrow = 3)
C[1,] = c("A", "B", "A", "A", "C")
C[2,] = c("B", "C", "A", "B", "C")
C[3,] = c("A", "B", "B", "A", "C")

seq = as.vector(C)
bpdf = prop.table(table(seq))
output = pssm(C, bpdf)

k = nrow(C)
sigma_size = length(unique(seq))
b = 1

#compute value for each count
q0 = (0+b)/(k + b*sigma_size)
q1 = (1+b)/(k + b*sigma_size)
q2 = (2+b)/(k + b*sigma_size)
q3 = (3+b)/(k + b*sigma_size)

q = as.matrix(cbind(c(q2,q1,q0), c(q0,q2,q1), c(q2,q1,q0), c(q2,q1,q0), c(q0,q0,q3) ))
p = matrix(c(0.4, 1/3, 4/15), ncol = 5, nrow = 3)
ans = q/p
rownames(ans) = c("A","B","C")
expect_equal(output, ans)

