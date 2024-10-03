#27th Feb 2024

#Script for turning every song into a spectogram, then excising the spectrograms for every note using
#the unit table from Lewis et. al. 2021. 

#path for all recordings
song_path = "~/Dropbox (The University of Manchester)/FINAL FILES/RenamedKoe/"
#get names of all song 
song_wavs = list.files(song_path, pattern = ".wav")
#grab only the wav files
song_wavs = song_wavs[grepl("\\.wav$", song_wavs)]

#read in song wav files individually and turn them into whole song spectrograms

song_specs = lapply(song_wavs, function(w){
  #for every song wav file, w do:
  
  #get path to the song
  path = paste(song_path, w, sep="")
  #read in wav file, starting from 0s
  w_wav = tuneR::readWave(path, units = "seconds", from = 0)
  #convert wav file to spectrogram, parameters come from Lewis et. al. 2021
  w_specgram = seewave::spectro(wave = w_wav, f = 44100, wl = 512, ovlp = 50, plot = FALSE) 
  
  #create output object
  output = w_specgram
  output$song_name = w
  return(output)
})

#need to extract specgrams for every note in unit table (from Lewis et. al. 2021)

#read in unit table
unit_tab <- read.csv( "~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv", header=TRUE )

#for every song recording (sound.files), we excise the spectrograms for every note
library(magrittr)

final_note_specgrams = lapply(song_wavs, function(w){
  #for every wav file, w do:
  song_tab = unit_tab %>%
    dplyr::filter(sound.files == w)
  
  #some songs haven't been processed and aren't in the unit table, since we can't use next in lapply
  if(nrow(song_tab)==0) return(NA)
  
  #get the whole song specgram
  song_specgram = song_specs[[which(song_wavs==w)]]
  
  #check correct specgram is loaded
  testthat::expect_equal(w, song_specgram$song_name)
  
  #cut spectrogram for every row in song_tab
  note_specgrams = vector(mode = "list", length = nrow(song_tab))
  for(r in 1:nrow(song_tab)){
    row = song_tab[r,]
    time = song_specgram$time
    #find the closest time cutoff, directly preceding note start
    cut_start = which.max(time[time<row$start])
    #find closest time cutoff, directly after note end
    cut_end = which(time>row$end)[1]
    
    note_spec = list(time = song_specgram$time[cut_start:cut_end], 
                     freq = song_specgram$freq,
                     amp = song_specgram$amp[,cut_start:cut_end],
                     song_name = w, 
                     selec = row$selec,
                     note_label = row$note_label,
                     song_individual = row$song_individual
                     )
    note_specgrams[[r]] = note_spec
  }
  
  return(note_specgrams)
})

final_note_specgrams = lapply(final_note_specgrams, function(x) x[!is.na(x)])

saveRDS(final_note_specgrams, file = "./data/final_note_specgrams.rds")

#compute sinkorn distances between note classes----

#small example
rec1 = final_note_specgrams[[1]]

note1 = rec1[[1]]
note2 = rec1[[2]]

#get amplitude matrices
A1 = note1$amp
A2 = note2$amp

#process amplitude matrices
As = lapply(list(A1,A2), function(X){
  #turn everything positive
  pos = X - min(X)
  #normalise
  pos/sum(pos)
})

#vectorise matrices for optimal transport ----
vA1 = as.vector(As[[1]])
vA2 = as.vector(As[[2]])

#Cost matrix
source("./functions/opt_functions/cost_mat.R")
C = cost_mat(p = note1, q = note2)

#compute sinkorn distances
setwd("./functions/opt_functions/")
source("SinkhornDistances.R")
res = entropyRegularisedKOT(p = vA1, q = vA2, cost.mat = C,max.cycles = 500)



####