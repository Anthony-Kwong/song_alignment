library("aphid")
data("globins")
globins

#phmm example

#fit phmm model
globins.PHMM <- derivePHMM(globins, residues = "AMINO", pseudocounts = "Laplace")
#plot phmm statespace
plot(globins.PHMM)
#use Viterbi to get path
path <- Viterbi(globins.PHMM, globins["GLB1_GLYDI", ])$path
#rename path in Delete, Match, Insert space
c("D", "M", "I")[path + 1]

#simulate sequences
sim <- list(length = 10)
suppressWarnings(RNGversion("3.5.0"))
set.seed(9999)
for(i in 1:10) sim[[i]] <- generate(globins.PHMM, size = 20)
sim

#we can optimise phmm 
globins2.PHMM <- train(globins.PHMM, sim, method = "BaumWelch", 
                       deltaLL = 0.01, seqweights = NULL)

#alignment using phmm
globins <- unalign(globins)
align(globins, model = globins.PHMM, seqweights = NULL, residues = "AMINO")
