#plot select alignments for results

library(ggplotify)

#make exemplar plots for every technique ----

#gibbs ----

gibbs_names = list.files("./results/fasta/robust_test/gibbs/")
lines = c("Blue", "Dark Pink", "Green", "Orange", "PaleBlue", "PaleGreen", "Pink", "Turquoise", "White", "Yellow")

for(i in 1:10){
  #setting indices
  start = 3*(i-1) + 1
  end = start + 2
  tur_gibbs = paste("./results/fasta/robust_test/gibbs/",gibbs_names[start:end], sep ="")
  gibb_plots = lapply(tur_gibbs, function(d){
    ggmsa(d, color="LETTER", seq_name=FALSE, char_width=0.4)
    #+ geom_msaBar() #insert consensus bar
  })
  gibb_plots = lapply(gibb_plots, as.ggplot)
  gibb_turq = ggpubr::ggarrange(plotlist = gibb_plots, ncol = 3)
  
  #set file name
  fname = paste("./results/alignment_plots/gibbs/gibbs_", lines[i],".png",sep = "")
  ggsave(plot = gibb_turq, file = fname)
}

#progressive alignment ----

dp_scores = readr::read_csv(file = "./results/msa_scores/dp_scores.csv")
pa_names = list.files("./results/fasta/robust_test/dynam_prog/")

#plot every alignment individually
for(i in seq_along(pa_names)){
  target = paste("./results/fasta/robust_test/dynam_prog/", pa_names[i], sep = "")
  
  #cleaned_target = read_fasta_sequences(target)
  
  p = ggmsa(target, color="LETTER", seq_name = FALSE, char_width=0.5)
  #+geom_msaBar()
  fname = paste("./results/alignment_plots/progress_align/",pa_names[i], ".png", sep ="")
  ggsave(plot = p, filename = fname)
}

#pHMM

phmm_scores = readr::read_csv("./results/msa_scores/phmm_scores.csv")
phmm_names = list.files("./results/fasta/robust_test/pHMM/")

for(i in seq_along(phmm_names)){
  target = paste("./results/fasta/robust_test/pHMM/", phmm_names[i], sep = "")
  p = ggmsa(target, color="LETTER", seq_name = FALSE, char_width=0.5)
  #+geom_msaBar()
  fname = paste("./results/alignment_plots/phmm/",phmm_names[i], ".png", sep ="")
  ggsave(plot = p, filename = fname)
}

#printing side by side alignments for comparison in chapter

# z = bird_songs %>% dplyr::filter(Line == "Pale Green")
# z$Bird.ID %>% unique()
# # #function to read fasta files, for debugging
# read_fasta_sequences <- function(path) {
#   raw <- readLines(path)
# 
#   # Remove header lines (starting with ">")
#   seqs <- raw[!startsWith(raw, ">")]
# 
#   # Remove whitespace and join into full strings if multi-line
#   seqs <- gsub("\\s+", "", seqs)
# 
#   # Return sequences as a pure character vector
#   return(seqs)
# }
# 
# path = paste("./results/fasta/robust_test/pHMM/", phmm_names[72], sep = "")
# z = read_fasta_sequences(path)

