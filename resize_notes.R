#We resample all the note spectrograms into the same size image (3d tensor)

note_specgram = readRDS(file = "~/work/sound_files/single_note_spectrograms/allnote_specgram.rds")
source("./functions/ResampleSpectrogram.R")

x = note_specgram[[1]]
y = resample.spectrogram(x, new.n.cols = 40)

#start 1558
#based on our eda, we resample everybody into 40 columns (approx mean of the number of cols)
res_specgram = lapply(note_specgram, function(x){resample.spectrogram(x, new.n.cols = 40)})

# z = res_specgram[[1]]
# image( t(z$amp), x=z$time, y=z$freq )
# amp = z$amp - min( z$amp)
# amp = amp/sum(amp)

#write loop to normalize the matrix so they are distributions
norm_specgram = lapply(res_specgram, function(x){
  amp = x$amp - min(x$amp)
  amp = amp/sum(amp)
  x$dist = amp
  return (x)
})

