# Given file name representing a note, extract fields from it usefully.
# Names look like JS0002-20110427-001-1.wav

library(stringr) 
library(chron)

parse.note.filename <- function( filename ) {
	# Split the name into dash-separated tokens
	token <- strsplit( filename, '-' )[[1]]
	if( length(token) != 4 ) {
		warning( sprintf( "Note file's name, %s, has unexpected format.", filename ))
		return( NA )
	}
	
	# Get the bird's number from the first token
	js.num <- strtoi( substr(token[1], 3, nchar(token[1])), base=10 )
	
	# Get the date from the second token
	year <- strtoi( substr(token[2], 1, 4), base=10 )
	month <- strtoi( substr(token[2], 5, 6), base=10 )
	day <- strtoi( substr(token[2], 7, 8), base=10 )
	date.str <- sprintf( "%d-%d-%d", day, month, year )
	recording.date <- chron( date.str, format='d-m-y' )
	
	
	# Get the song number from the third token
	song.num <- strtoi( token[3], base=10 )
	
	# Finally, the note number from the 4th token
	note.num <- strtoi( str_extract(token[4], "^\\d+"), base=10 )
	
	# Assemble the result as a named list and return
	return( list(ID=js.num, date=recording.date, song=song.num, note=note.num) )
	
}



name.list <- c( "JS0002-20110427-001-1.wav",  
	"JS0002-20110427-001-10.wav", "JS0002-20110427-001-11.wav",
	"JS0002-20110427-001-12.wav", "JS0002-20110427-001-13.wav",
	"JS0002-20110427-001-14.wav", "JS0002-20110427-001-15.wav",
	"JS0002-20110427-001-16.wav", "JS0002-20110427-001-17.wav",
	"JS0002-20110427-001-2.wav"  
)

parse.note.filename( name.list[1] )
sapply( name.list, parse.note.filename )
