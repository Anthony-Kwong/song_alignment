#We resample all the note spectrograms into the same size image (3d tensor)

note_specgram = readRDS(file = "~/work/sound_files/single_note_spectrograms/allnote_specgram.rds")
source("./functions/ResampleSpectrogram.R")

x = note_specgram[[1]]
y = resample.spectrogram(x, new.n.cols = 40)

#start 1558
#based on our eda, we resample everybody into 40 columns (approx mean of the number of cols)
res_specgram = lapply(note_specgram, function(x){resample.spectrogram(x, new.n.cols = 40)})
