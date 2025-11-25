#generate a pairwise alignment of a pupil and a tutor as an illustrative example

library(aphid)
library(ggmsa)
library(magrittr)
pacman::p_load(aphid, ggmsa, magrittr, ggplot2)

#load function
source("./functions/generate_alignment.R")

bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")

#align 2 birds
sub_birds = bird_songs %>%
  dplyr::filter(bird.num %in% c(37, 52))

b_fasta = generate_alignment(sub_birds)
fname = "./seq_evo/paper/results/pupil_tutor.fasta"
bio3d::write.fasta(b_fasta, file = fname)

plot = ggmsa(fname , color = "LETTER", seq_name = TRUE, char_width = 0.2) + geom_msaBar() 

plotname = paste("./seq_evo/paper/results/pupil_tutor.pdf")
ggsave(plot, file = plotname)
