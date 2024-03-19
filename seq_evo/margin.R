#script to compare how much transport is happening on the frequency dimension vs time
#since there's many more frequency bins relative to time bins.

#read in spectrograms, every element of list is a song
final_note_specgrams <- readRDS("./data/final_note_specgrams.rds")
#un-nest such that element list element is now a note
final_note_specgrams <- unlist(final_note_specgrams, recursive = FALSE)

#
# x = final_note_specgrams[[1]]
# x[2]

#create get functions for ease of value retrieval 

#get_vals ----

#input: specgram object
#value : the note label and song individual
get_vals <- function(specgram){
  #x = unlist(specgram, recursive = F)
  output = tibble::tibble(note_label = x$note_label, song_individual = x$song_individual)
  return(output)
}

#small example ----

# #select birds
# b1 = "JS0002"
# b2 = "JS0003"
# 
# #retrieve note labels and singer for every specgram
# k = lapply(final_note_specgrams, get_vals)
# specgram_tab = do.call(rbind, k)
# 
# #get most common note label
# b_tab = specgram_tab %>%
#   dplyr::filter(song_individual == b1)
# table(b_tab)
# 
# #get indicies of the b1
# b_index = which(specgram_tab$song_individual == b1)
# 
# n1 = final_note_specgrams[[1]]
# n2 = final_note_specgrams[[2]]
# 
# source("./functions/opt_functions/margin_opt.R")
# z = margin_opt(n1,n2)

#compare random samples of note pairs
source("./functions/opt_functions/margin_opt.R")

nsam = 10000
S = seq(length(final_note_specgrams))
n1s = sample(S, nsam, replace = T)
n2s = sample(S, nsam, replace = T)

fres = rep(NA, nsam)
tres = rep(NA, nsam)

for(i in 1:nsam){
  n1 = final_note_specgrams[[n1s[i]]]
  n2 = final_note_specgrams[[n2s[i]]]
  test = margin_opt(n1,n2)
  fres[i] = test$freq$optimal.cost
  tres[i] = test$time$optimal.cost
}

x = tibble::tibble(cost = fres, dim = "freq")
y = tibble::tibble(cost = tres, dim = "time")
res_tab = rbind(x,y)

median(fres)/median(tres)

library(ggplot2)
ggplot(res_tab, aes(x = dim, y = cost)) + 
  geom_boxplot()

coef_f = sd(fres)/mean(fres)
coef_t = sd(tres)/mean(tres)
mean(fres)
mean(tres)

readr::write_csv(res_tab, file = "./data/margin.csv")
