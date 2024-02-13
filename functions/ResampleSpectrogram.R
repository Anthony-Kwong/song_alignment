source( "./functions/resampleStepFunc.R")

#resample spectrogram function

#sg: current spectrogram (list of numeric vectors with time,freq,amp)
#new.n.cols : Numeric scalar for number of new boundary points 

resample.spectrogram <- function( sg, new.n.cols )
{
	orig.n.cols <- length(sg$time)
	orig.bdy <- (0:orig.n.cols) / orig.n.cols
	new.bdy <- (0:new.n.cols) / new.n.cols
	new.amp <- apply( 
		sg$amp, MARGIN=1, # by rows
		FUN = function( row ) {
			orig.sf <- list( bdy=orig.bdy, val=row )
			new.sf <- resample.step.func( orig.sf, new.bdy )
			return( new.sf$val )
		}
	)
	
	new.sg <- sg
	new.sg$amp <- t(new.amp)
	new.sg$time <- seq( from=sg$time[1], to=sg$time[orig.n.cols], length.out=new.n.cols)
	return( new.sg )
}

#examples
# note_specgram = readRDS(file = "~/work/sound_files/single_note_spectrograms/allnote_specgram.rds")
# allnote_specgram = note_specgram[[1]]
# new.sg <- resample.spectrogram( allnote_specgram, 5)
# image( t(allnote_specgram$amp), x=allnote_specgram$time, y=allnote_specgram$freq )
# image( t(new.sg$amp), x=new.sg$time, y=new.sg$freq )


