#exploratory analysis on note sequences

song.df <- read.csv( "~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv", header=TRUE )
summary( song.df )

# Drop everything but file and note
small.df <- subset( song.df, TRUE, select=c("sound.files","note_label", "song_individual") )
names(small.df ) <- c( "sound.file","note.label", "song_individual" )
small.df$sound.file <- factor(small.df$sound.file)
head(small.df )

################################################################
#   Give the notes single-character names, chosen so that
# the most common ones get letters earlier in the alphabet
################################################################

# Count the number of occurrences of each note label
count.df <- xtabs( ~song_individual + note_label, data = song.df)
count.df

stats::heatmap(t(count.df))
gplots::heatmap.2(t(count.df), col = terrain.colors(256))


