#first step in likelihood work

library(magrittr)

bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
meta.data = readr::read_csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Files for Anthony/JavaSparrow_Metadata.csv")

#load functions
source("./functions/align_mat.R")
source("./functions/norm_pairdiff.R")
source("./functions/dist_mat.R")

#example with pink lineage
k = bird_songs %>%
  dplyr::filter(Line == "Blue")
birds = unique(k$Bird.ID)

#align bird pairs
bird_comb = as.data.frame(t(combn(birds, m = 2)))
pairs = nrow(bird_comb)

#intialise vector of distances
D = rep(NA, pairs)

for(i in 1:pairs){
  p = bird_comb[i,]
  #align songs
  ap = align_mat(p$V1, p$V2)
  #extract alignment matrices
  As = ap$alignment
  #make pairwise comparisons
  D[i] = norm_pairdiff(s1 = As[[1]],As[[2]])
}

#add the distances to to main dataframe
res = cbind(bird_comb, D)
bird_dist = dist_mat(res)
tre = nj(bird_dist)
plot(tre, "u")

#Method 1: Adapt the "morphological" method. We used normalized pairwise differences.

### From Saitou and Nei (1987, Table 1):
library(ape)
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(0, 8, 8)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(1:8, 1:8)
tr <- nj(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- nj(dist.dna(woodmouse))
plot(trw)


