library( stringr )
library( collections )

# Read the note table
song.df <- read.csv( "JavaSparrow_UnitTable.csv", header=TRUE )
summary( song.df )

# Drop everything but file and note
small.df <- subset( song.df, TRUE, select=c("sound.files","note_label") )
head( small.df )

# Give the notes shorter names: we do this with a dictionary (a.k.a. a hashmap)
n.notes <- length( unique(small.df$note_label) )
note.to.char <- dict( keys=unique(small.df$note_label), items=LETTERS[1:n.notes] )
with( note.to.char, cbind(keys(), values() )) # sanity check

small.df$note.char <- sapply( small.df$note_label, function(str) { note.to.char$get(str) })
head( small.df )

# Concatenate all the notes in a given song into a single string
seq.df <- aggregate( note.char ~ sound.files, 
                     FUN=function(c) { str_c(c, sep="", collapse="") }, 
                     data=small.df)

names( seq.df ) <- c( "sound.files", "note.seq")
head( seq.df )
