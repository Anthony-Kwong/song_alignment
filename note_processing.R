#script to generate single wav files and single spectrograms for all notes
library("warbleR")

setwd("~/Dropbox (The University of Manchester)/FINAL FILES/RenamedKoe/")
A = list.files()
B = list.files(pattern = ".jpeg")

#get names of wav files only
wavs = A[!A %in% B]

#read in unit table
unit_tab <- read.csv( "~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv", header=TRUE )

#we filter out notes that are too short because we don't want to pad them, it's about 2percent of chirps

#44kHz samples so...... 44*unit_tab$duration>512, i.e. duration cutoff is 512/44

unit_tab = unit_tab %>%
  dplyr::filter(duration > 512/44)

#cut selections (notes) into individual sound file, using the unit table and the full recordings in Dropbox----

#we added padding (mar) because some notes are really short
cut_sels(X = unit_tab, mar = 0 , path = "~/Dropbox (The University of Manchester)/FINAL FILES/RenamedKoe/", 
         dest.path = "~/work/sound_files/single_note_wavs/")

#get names of the new wavs
note_wavs = list.files("~/work/sound_files/single_note_wavs/")
#check the counts agree
nrow(unit_tab) == length(note_wavs)
#check names agree
unit_tab$sound.files[5]
note_wavs[5]

#turn wavs into spectrograms, check the options Becky used in Frontiers paper when she used KoE
# setwd("~/work/sound_files/single_note_wavs/")
# eg = "JS0002-20110427-001-1.wav"
# file = read_wave(X = eg)
# output = seewave::spectro(wave = file, wl = 512, ovlp = 50)

source("~/Documents/GitHub/song_alignment/functions/split_string.R")

problem = "JS0003-20110501-011-42.wav"

#change to note directory
setwd("~/work/sound_files/single_note_wavs/")

#convert the wav files in note_wavs into spectrograms via seewave::spectro. Add in metadata via unit_tab
spectrograms = lapply(note_wavs, function(w){
  print(w)
  #add meta data 
  note_identity = split_string(w)
  #find row corresponding to the sound.files and selec
  
  #get mini unit table of the recording based on sound.files
  song_tab = unit_tab[which(unit_tab$sound.files == note_identity$sound.files),]
  
  #get the correct note using selec
  note_row = song_tab[which(song_tab$selec == note_identity$selec),]
  
  #read in note as wav file
  file = read_wave(w)
  #parameters for spectro were taken from Becky Frontiers paper, still need to check dB
  #dB weights https://en.wikipedia.org/wiki/A-weighting, https://en.wikipedia.org/wiki/ITU-R_468_noise_weighting
  #we used default after checking KoE code
  
  output = seewave::spectro(wave = file, wl = 512, ovlp = 50, plot = FALSE) #get everything on a linear scale
  output$sound.files = note_identity$sound.files
  output$selec = note_identity$selec
  output$note_label = note_row$note_label
  output
})

saveRDS(spectrograms, file = "~/work/sound_files/single_note_spectrograms/allnote_specgram.rds")


#we need to turn spectrograms into 2dim densities, notes have different durations----

#read in spectrograms
note_specgram = readRDS(file = "~/work/sound_files/single_note_spectrograms/allnote_specgram.rds")
x = note_specgram[[1]]

#sanity check for spectrograms
eg = "JS0002-20110427-001-10.wav"
file = read_wave(X = eg)
output = seewave::spectro(wave = file, wl = 512, ovlp = 50)

#every specgram is an object 


