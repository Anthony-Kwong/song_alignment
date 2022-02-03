################################################################
# Preprocessing to enable a first stab at aligning songs
# 
# mrm and Anthony Kwong: Man Uni, 9 December 2021
################################################################
rm( list=ls() ) # wipe the slate clean

library( stringr ) # str_c()
library( collections ) # dict()

################################################################
# Read Becky's unit table, then drop unnecessary columns
################################################################
song.df <- read.csv( "~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv", header=TRUE )
summary( song.df )

# Drop everything but file and note
small.df <- subset( song.df, TRUE, select=c("sound.files","note_label") )
names(small.df ) <- c( "sound.file","note.label" )
head(small.df )

################################################################
#   Give the notes single-character names, chosen so that
# the most common ones get letters earlier in the alphabet
################################################################

# Count the number of occurrences of each note label
count.df <- aggregate( sound.file ~ note.label, length, data=small.df )
names( count.df ) <- c( "note.label", "n.occurrences" )

# Get an ordered list of note labels
perm <- order( count.df$n.occurrences, decreasing=TRUE )
ordered.note.labels <- count.df$note.label[perm]

# Give the notes shorter names: we do this with a dictionary (a.k.a. a hashmap)
n.notes <- length( unique(small.df$note.label) )
note.to.char <- dict( keys=ordered.note.labels, items=LETTERS[1:n.notes] )

# Check that all has turned out as we intended
count.df <- count.df[perm,]
count.df$note.char <- sapply(count.df$note.label, 
                             function(str) { note.to.char$get(str) } 
)

count.df

# Tranlate note names in the main data frame
small.df$note.char <- sapply( small.df$note.label, 
                              function(str) { note.to.char$get(str) }
)

head( small.df )

################################################################
#   Finally, concatenate all the notes from a given
# recording into a single string.
################################################################

# Concatenate all the notes in a given song into a single string

#aggregate splits data by sound.file, then takes all the note.char and collapse
#them into a string
seq.df <- aggregate( note.char ~ sound.file, 
                     FUN=function(c) { str_c(c, sep="", collapse="") }, 
                     data=small.df)

names( seq.df ) <- c( "sound.file", "note.seq")
head( seq.df )

################################################################
#   Parse the file names to extract dates as well as
# bird and recording numbers
################################################################

# The names are of the form "JS0002-20110427-001.wav"
# We being by splitting them into dash-separated tokens
fnameTokenList = sapply( seq.df$sound.file, 
                         function(fname) { strsplit(fname, '-'  ) }  
)

# Get the bird's number
birdNums <- vapply( fnameTokenList, 
                    function(tokens){
                      birdStr <- tokens[1]
                      birdNumStr <- gsub( 'JS[0]+', '', birdStr, perl=TRUE )
                      return( as.integer(birdNumStr) )
                    }, 
                    FUN.VALUE=c(birdNum=0), # Specify the return type
                    USE.NAMES=FALSE
)

# Get a string representing the date
recordingDates <- vapply( fnameTokenList, 
                          function(tokens){
                            rawDateStr <- tokens[2]
                            year <- substr( rawDateStr, 1, 4 )
                            month <- substr( rawDateStr, 5, 6 )
                            day <- substr( rawDateStr, 7, 8 )
                            dateStr <- paste( day, month, year, sep="/" )
                            return( dateStr )
                          }, 
                          FUN.VALUE=c(date="01/01/2011"), 
                          USE.NAMES=FALSE
)

# Extract the recording number
recordingNums <- vapply( fnameTokenList, 
                         function(tokens){
                           myStr <- tokens[3]
                           recNumStr <- gsub( '.wav', '', myStr, perl=TRUE )
                           return( as.integer(recNumStr) )
                         }, 
                         FUN.VALUE=c(recNum=0), 
                         USE.NAMES=FALSE
)

# Add the columns
seq.df$bird.num <- birdNums
seq.df$rec.num <- recordingNums
seq.df$rec.date <- recordingDates
head( seq.df )

#Now we are ready to do some alignments

#try some distance matrices




