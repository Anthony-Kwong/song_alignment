#script to investigate consensus scores in alignments

#load functions
source("./functions/fasta2matrix.R")
source("./functions/consensus_score.R")
source("./functions/shannon_score.R")

#small example
# ex = fasta2matrix("./data/songs_fasta/JS100.fasta")
# A = get_column_proportions(ex)
# res = consensus_score(A)

#compute consensus for every lineage
setwd("./results/fasta/lines_fasta/")
lins = list.files()

res = lapply(lins, function(d){
  #read in fasta 
  A_mat = fasta2matrix(d)
  prop_mat = get_column_proportions(A_mat)
  #compute consensus
  con = consensus_score(prop_mat)
  #compute shannon
  shannon = shannon_score(prop_mat)
  #change the name to remove .fasta
  name = sub("\\.fasta$", "", d)
  tibble::tibble(Lineage = name, consensus = con, shannon = shannon)
})
res = do.call(rbind, res)

library(ggplot2)

#change colours to reflect those in the name
colours <- c("Blue" = "blue", "Dark Pink" = "deeppink", "Green" = "green", 
             "Orange" = "orange", "Pale Blue" = "lightblue", 
             "Pale Green" = "lightgreen", "Pink" = "pink", 
             "Turquoise" = "turquoise", "White" = "white", "Yellow" = "yellow")

con_plot = ggplot(data = res, aes(y = consensus, fill = Lineage)) +
  geom_boxplot() +
  scale_fill_manual(values = colours) +
  theme(
    legend.position = "bottom",                         # Move legend below
    legend.direction = "horizontal",                    # Arrange legend horizontally
    axis.text.x = element_blank(),       
    axis.ticks.x = element_blank()
  )

ggsave(plot = con_plot, filename  = "../../../results/consensus_plot.png")

#probably won't use shannon

ggplot(data = res, aes(y = shannon, fill = Lineage)) +
  geom_boxplot() +
  scale_fill_manual(values = colours) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom",                         # Move legend below
    legend.direction = "horizontal"                     # Arrange legend horizontally
  ) + 
  ylab("shannon entropy")


