################################################################
# This is the master preprocessing script which turns the original unit table Lewis. et.al 
# into note sequences (NoteSequences.csv). 
# We also generate a heatmap of note usages. 
# 
# mrm and Anthony Kwong: Man Uni, 9 December 2021
################################################################
rm( list=ls() ) # wipe the slate clean

library( stringr ) # str_c()
library( collections ) # dict()
library(birdsong.tools)

################################################################
# Read Becky's unit table, then drop unnecessary columns
################################################################
song.df <- read.csv( "~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv", header=TRUE )
summary( song.df )
song.df <- dplyr::rename(song.df, Bird.ID = song_individual)
meta.data=read.csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Files for Anthony/JavaSparrow_Metadata.csv")
#issue with add_metadata (can't tolerate single column)
song.df <- birdsong.tools::add_metadata(song.df, meta.data, cols = c(5,6))

data = tibble::tibble(Bird.ID = c("JS001", "JS002"),a = c(1,2), b = c(1,2))
metadata2 = tibble::tibble(Bird.ID = c("JS001", "JS002"), x= c(3,4), y = c(5,6), z = c(7,8))
add_metadata(data,metadata2, cols = 3)

# Drop everything but file and note
small.df <- subset( song.df, TRUE, select=c("sound.files","note_label","Clutch", "Bird.ID") )
names(small.df ) <- c( "sound.file","note.label","clutch","Bird.ID")
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

#plot bar chart of general note usage
library(ggplot2)
ggplot2::ggplot(data = count.df, aes(y = n.occurrences, x = note.label)) +
  geom_bar(stat = "identity")

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
seq.df <- aggregate( note.char ~ sound.file + clutch, 
                     FUN=function(c) { str_c(c, sep="", collapse="") }, 
                     data=small.df)

names( seq.df ) <- c( "sound.file","clutch", "note.seq")
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

singer = vapply( fnameTokenList, 
                 function(tokens){
                   birdStr <- tokens[1]
                   return( birdStr )
                 }, 
                 FUN.VALUE=c(birdNum="string"), # Specify the return type
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
seq.df$Bird.ID <- singer
head( seq.df )

#insert line info
seq.df = add_metadata(seq.df, meta.data, cols = c(5:9))

################################################################
#   Blap everything out to files
################################################################

write.csv( count.df, "./data/NoteNames.csv", row.names=FALSE )
write.csv( seq.df, "./data/NoteSequences.csv", row.names=FALSE )

################################################################
#   Do some EDA ----
################################################################
library(gplots)
# Look at note usage in the sense of raw note count
note.count.mat <- xtabs( ~ Bird.ID + note_label, data=song.df )
pdf( "NoteUsageCountHeatmap.pdf", height=210/25.4, width=297/25.4 )
heatmap.2( 
  t(note.count.mat), 
  col=terrain.colors(256), trace="none",
  main="Number of notes recorded",
  key.title="Key", key.xlab="Notes recorded",
  margins=c(4,7)
)
dev.off()

#Turn counts into proportions because some birds might just sing more notes. 

# Compute the fraction of each bird,s notes that are of each type
notes.per.bird <- rowSums( note.count.mat )
#compute proportions
note.prop.mat <- diag( 1.0 / notes.per.bird ) %*% note.count.mat 
rownames( note.prop.mat ) <- rownames( note.count.mat ) # Names seem to get lost
rowSums(note.prop.mat ) # Should all be 1.0

# Plot new heatmap based on proportions
pdf( "NoteUsageProportionHeatmap.pdf", height=210/25.4, width=297/25.4 )
heatmap.2( 
  t(note.prop.mat), 
  col=terrain.colors(256), trace="none",
  main="Proportion of notes recorded (per bird)",
  key.title="Key", key.xlab="Proportion of notes",
  margins=c(4,7)
)
dev.off()

#generate colors for each clutch variable
clutches = unique(seq.df$clutch)
clutch_cols = rainbow(length(clutches))
clutchd <- dict(clutch_cols, clutches)

clutchd$get("A")

x_colors = sapply(rownames(note.prop.mat), function(x){
  #find the clutch of each individual
  clutch = metadata[which(metadata$Bird.ID == x),]$Clutch
  #find the color for the given clutch
  clutchd$get(clutch)
})

pdf( "NoteUsageProportionHeatmap.pdf", height=210/25.4, width=297/25.4 )
heatmap.2( 
  t(note.prop.mat), 
  col=terrain.colors(256), trace="none",
  main="Proportion of notes recorded (per bird)",
  key.title="Key", key.xlab="Proportion of notes",
  margins=c(4,7), 
  colCol = x_colors
)
dev.off()

#number of notes
seq.df[1,]$note.seq %>% nchar()

seq.df <- seq.df %>% 
  dplyr::mutate(num_notes = nchar(note.seq))

#histogram of note number per sequence
ggplot2::ggplot(seq.df, aes(num_notes)) + 
  geom_histogram(bins = 10)

#ggplot heatmap example

m <- matrix(rnorm(400), ncol=40)
sample.types <- c(rep("Blue", 10), rep("Green", 10), rep("Red", 10), rep("Purple", 10))

library(gplots)
heatmap.2(m, trace="none", colCol = sample.types)

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)
data$class <- sample(c("A","B"),20,replace = T)

# Heatmap 
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()

## adding new abbreviated note labels to the unit table for convenience
head(song.df)

#couldn't get dict to work in mutate
# song.df %>%
#   dplyr::mutate(label2 = note.to.char$get(note_label))

#work around
new_names = lapply(song.df$note_label, function(x){note.to.char$get(x)}) %>%
  unlist()

song.df2 = tibble::tibble(song.df, notelab2 = new_names)
write.csv(song.df2, file = "./data/Unit_tab_dict.csv")

##Alignments ----

#Compute pairwise distances between sequences

# library(alakazam)
# 
# pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
#              #need an input distance matrix to show closeness of different individual elements
#              dist_mat=getDNAMatrix(gap=0))
# 
# #Try out some classic alignment algorithms (e.g. Needleman Wunsch)
# 
# #check by hand
# 
# head(seq.df)
# dist(seq.df)
# 
# #Now we are ready to do some alignments
# 
# #examples seqs
# ex_songs = seq.df$note.seq[1:5]
# pairwiseDist(ex_songs)
# 
# ##edit distances
# library(stringdist)
# stringdist(ex_songs[1], ex_songs[2], method = 'lv')
# 
# 
# d <- stringdistmatrix(c('foo','bar','boo','baz'))
# # try 
# plot(hclust(d))
# 
# stringdistmatrix(c("foo","bar","boo"),c("foo","bar","boo"))


#need to add gaps in before we can compute differences

#try some distance matrices then alignment algorithms

#steps to do

#1. Compute edit distances (don't need alignment I think)
#2. attempt to align character strings using basic alignment algorithnms
#3. Do pairwise distances based on alignments
#4. Substitution models and compute likelihoods of trees


#Alignment helps us to identify point mutations and indels
