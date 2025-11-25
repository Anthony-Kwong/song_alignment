#Script for classifying songs using pHMMs

#load libs
library(magrittr)
library(aphid)
library(reactable)

#load functions
source("./functions/string_split.R")
source("./functions/bck_pdf.R")

#read in whole song data set
bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
lines = unique(bird_songs$Line)

#store classification results 
line_test_res = list()

#loop thru lineages
for(l in 1:length(lines)){
  
  #filter down dataset to lineage
  line_data = bird_songs %>%
    dplyr::filter(Line == lines[l])
  
  #retrieve alphabet ----
  
  #get note sequences as long strings
  bird_songsseqs = line_data$note.seq
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  #combine all strings together
  comb = paste(bird_songsseqs, collapse = "")
  #take all unique letters
  alphabet = unique(strsplit(comb, "")[[1]])
  
  #preallocate for testing data
  line_test = list()
  #store model for every singer
  line_mods = list()
  
  singers = unique(line_data$Bird.ID)
  #loop for every singer
    for(s in 1:length(singers)){
      singer_data = line_data %>%
        dplyr::filter(Bird.ID == singers[s]) %>%
        dplyr::select(Bird.ID,note.seq, Line)
      
      nsongs = nrow(singer_data)
      #partition into training and testing set
      sam = sample(seq(nsongs), size = ceiling(nsongs/2))
      train = singer_data[sam,]
      line_test[[s]] = singer_data[-sam,]
      
      #model training ----
      
      #processing note sequences into split strings
      train_seqs =string_split(train$note.seq)
      #train pHMM 
      line_mods[[s]] <- derivePHMM(train_seqs, residues = alphabet, pseudocounts = "Laplace", refine = "BaumWelch")
    }
  #classification for singers in lineage ----
  
  #add names to the models by singer
  names(line_mods) = singers
  
  
  #bind all testing data together
  final_test = do.call(rbind, line_test)
  #store predictions
  bird_pred = rep(NA,nrow(final_test))
  #for every row, put note.seq through every model
  for(r in 1:nrow(final_test)){
    bird_row = final_test[r,]
    #get song for testing
    test_song = string_split(bird_row$note.seq)
    #loop through every model pHMM
    
    bird_scores = rep(NA, length(line_mods))
    for(m in 1:length(line_mods)){
      #get Viterbi paths and grab scores for every model
      path = Viterbi(line_mods[[m]], test_song)
      bird_scores[m] = path$score
    }
    print(bird_scores)
    #make prediction on row r based on highest path score
    best_scorer = which.max(bird_scores)
    bird_pred[r] = names(line_mods)[best_scorer]
  }
  #save lineage test results
  line_test_res[[l]] = cbind(final_test, prediction = bird_pred)
}

#bind all results together and write csv
singer_class_res = do.call(rbind, line_test_res)
#add T/F column
singer_class_res = singer_class_res %>%
  dplyr::mutate(correct = ifelse(prediction == Bird.ID, TRUE, FALSE))
#compute missclass rate (whole population)
sum(singer_class_res$correct)/nrow(singer_class_res)
#save
readr::write_csv(singer_class_res, file = "./results/song_class/song_class.csv")

#compute missclass rate for every lineage separately

class_tab = singer_class_res %>%
  dplyr::group_by(Line) %>%
  summarise(
    #fix nrow bit
    accuracy = sum(correct)/n(),
    N = n(),
    Birds = length(unique(Bird.ID)),
    base_rate = 1/length(unique(Bird.ID))
  )

xtable::xtable(class_tab)

#xtable::xtable(x, caption = "Table of note classes with their abbreviated name and number of occurences in the Java sparrows song dataset",                label = "note_tab", digits = 0)
