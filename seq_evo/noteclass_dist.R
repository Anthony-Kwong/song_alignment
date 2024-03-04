final_note_specgrams = readRDS(final_note_specgrams, file = "./data/final_note_specgrams.rds")

#un-nest such that element list element is now a note
final_note_specgrams <- unlist(final_note_specgrams, recursive = FALSE)

#compute Sinkorn distances between note classes----

#set sampling params
nsam = 3
nnotes = length(final_note_specgrams)

#obtain sample pairs
ps = sample(nnotes, size = nsam, replace = T)
qs = sample(nnotes, size = nsam, replace = T)

#source function
setwd("./functions/opt_functions/")
source("./note_compare.R")

res <- vector(mode = "list", length = nsam)

for(i in 1:nsam){
  n1 = final_note_specgrams[[ps[i]]]
  n2 = final_note_specgrams[[qs[i]]]
  
  #set vector to record convergence
  conv = rep(F, nsam)
  
  opt = note_compare(n1,n2)
  #record convergence result
  if(opt$converged == FALSE){
    conv[i] = T
  }
  
  #record note data for traceability
  opt$p_id = ps[i]
  opt$q_id = qs[i]
  opt$p_class = n1$note_label
  opt$q_class = n2$note_label
  
  res[[i]] = opt
}
