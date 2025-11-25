#script for generating note usage data for alignment paper

seq.df = readr::read_csv("./data/NoteSequences.csv")
#read in unit table with updated note labels
song.df <- read.csv( "./data/Unit_tab_dict.csv", header=TRUE )

#count note usage for every bird
note.count.mat <- xtabs( ~ Bird.ID + notelab2, data=song.df )

#Turn counts into proportions because some birds might just sing more notes. ----
# Compute the fraction of each bird,s notes that are of each type
notes.per.bird <- rowSums( note.count.mat )
#compute proportions
note.prop.mat <- diag( 1.0 / notes.per.bird ) %*% note.count.mat 
rownames( note.prop.mat ) <- rownames( note.count.mat ) # Names seem to get lost
rowSums(note.prop.mat ) # Should all be 1.0

#produce heatmap of note usage, where the cells are colored by the note proportion, birds
#are colored by their social lineage
# library(gplots)
# heatmap.2( 
#   t(note.prop.mat), 
#   col=terrain.colors(256), trace="none",
#   #main="Proportion of notes recorded (per bird)",
#   key.title="Key", key.xlab="Proportion of notes",
#   margins=c(4,7), 
#   colCol = x_colors,
#   dendrogram = "none"
# )

#want to put social lineages together....

meta.data=read.csv("~/Dropbox (The University of Manchester)/FINAL FILES/20210303/Files for Anthony/JavaSparrow_Metadata.csv")
rownames(note.count.mat)
meta.data$Bird.ID

library(dplyr)

birds = meta.data %>%
  dplyr::select(Bird.ID, Line)

birds_reordered = birds %>%
  arrange(Line, Bird.ID)

#reorder note proportions matrix rows
note.prop.mat = note.prop.mat[birds_reordered$Bird.ID, , drop = FALSE]

#generate colors for every social lineage----
sls = unique(seq.df$Line)
#generate colors for the lines, we can make this correspond to the names directly afterwards
sls_cols = colors <- c("blue", "deeppink", "turquoise", "black", "darkgoldenrod", "lightblue", 
                       "green", "hotpink", "orange", "palegreen")
slsd <- collections::dict(sls_cols, sls)

#example of getting value using key
slsd$get("Blue")

#use apply to get the colors for every bird in the note proportions matrix
x_colors = sapply(rownames(note.prop.mat), function(x){
  #find a row of every individual x
  ind_row = seq.df[which(seq.df$Bird.ID == x),][1,]
  line = ind_row$Line
  #find the color for the given clutch
  slsd$get(line)
})

library(gplots)

# Calculate the margins to center the heatmap
# heatmap_width <- ncol(note.prop.mat) * 0.2  # Assuming each column is 0.2 inches wide
# heatmap_height <- nrow(note.prop.mat) * 0.2  # Assuming each row is 0.2 inches tall
# left_margin <- (15 - heatmap_width) / 2
# bottom_margin <- (15 - heatmap_height) / 2

lwid=c(0.2,5) #make column of dendrogram and key very small and other colum very big 
lhei=c(0.2,5)

png("./results/eda/heatmap.png", width = 30, height = 30, units = "in", res = 300)
heatmap.2(
  t(note.prop.mat),
  col=terrain.colors(256), 
  trace="none",
  #main="Proportion of notes recorded (per bird)",
  key.title="Key", key.xlab="Proportion of notes",
  #reduce white margins
  margins=c(20, 20),
  colCol = x_colors,
  cexCol = 3, #column text sizes
  cexRow = 3,
  dendrogram = "none",
  #get R to not reorder rows and cols by dendrogram
  Rowv = NA, 
  Colv = NA, 
  key = F,
  #get rid of white spaces
  lwid = lwid, 
  lhei = lhei
)
dev.off()

blue_red_palette <- colorRampPalette(c("blue", "white", "red"))
#heatmap v2 with the color key
png("./results/eda/heatmap2.png", width = 45, height = 40, units = "in", res = 300)
heatmap.2(
  t(note.prop.mat),
  #col=terrain.colors(256), 
  col = blue_red_palette(256),
  trace="none",
  #main="Proportion of notes recorded (per bird)",
  key.title="Key", key.xlab="Proportion of notes",
  #reduce white margins
  margins=c(20, 20),
  colCol = x_colors,
  cexCol = 4, #column text sizes
  cexRow = 4,
  dendrogram = "none",
  #get R to not reorder rows and cols by dendrogram
  Rowv = NA, 
  Colv = NA, 
  key = F,
  #get rid of white spaces
  lwid = lwid, 
  lhei = lhei
)
dev.off()

#version 3 with scale bar

png("./results/eda/heatmap3.png", width = 12, height = 12, units = "in")
heatmap.2(
  t(note.prop.mat),
  #col=terrain.colors(256), 
  col = blue_red_palette(256),
  trace="none",
  #main="Proportion of notes recorded (per bird)",
  key.title="Key", key.xlab="Proportion of notes",
  #reduce white margins
  #margins=c(20, 20),
  colCol = x_colors,
  cexCol = 6, #column text sizes
  cexRow = 6,
  dendrogram = "none",
  #get R to not reorder rows and cols by dendrogram
  Rowv = NA, 
  Colv = NA, 
  key = T, #scale bar
  #keysize = 5
  #get rid of white spaces
  # lwid = lwid,
  # lhei = lhei
)
dev.off()


#boxplot of sequence length

lines = unique(seq.df$Line)
len_tab = lapply(lines, function(L){
  #filter for line
  line_df = seq.df %>% 
    dplyr::filter(Line == L)
  songs = line_df$note.seq
  #get song lengths by number of notes
  song_len = nchar(songs)
  #return dataframe
  tibble::tibble(note_count = song_len, Line = L)
})

len_tab = do.call(rbind, len_tab)

library(ggplot2)

dur_plot = ggplot(len_tab, aes(x = Line, y = note_count)) +
  geom_boxplot() +
  labs(x = "Social lineage", y = "song length")

ggsave(dur_plot, file = "./results/eda/songlen_boxplot.png")

#produce dictionary table
note_dic = readr::read_csv("./data/NoteNames.csv")
x = tibble::tibble(note_dic$note.label, note_dic$note.char, note_dic$n.occurrences)
colnames(x) = c("note name", "letter", "number of occurences")
xtable::xtable(x, caption = "Table of note classes with their abbreviated name and number of occurences in the Java sparrows song dataset",
               label = "note_tab", digits = 0)

#number of unique note classes per social lineage
head(seq.df)
lines

x = seq.df %>%
  dplyr::filter(Line == L)
y <- paste0(x$note.seq,sep = "", collapse = "")
strsplit(y, "")[[1]]
strsplit(y, "")[[1]] %>% unique()

char_tab = lapply(lines, function(L){
  ldf = seq.df %>%
    dplyr::filter(Line == L)
  lchar = paste0(ldf$note.seq,sep = "", collapse = "")
  #get number of unique characters
  nchar = strsplit(lchar, "")[[1]] %>% unique()
  unchar = length(nchar)
  
  #get number of songs
  nsong = nrow(ldf)
  
  #get number of birds
  nbird = unique(ldf$Bird.ID) %>% length()
  
  #output data
  tibble::tibble(Line = L, birds = nbird, nsong ,unchar)
})

char_tab = do.call(rbind, char_tab)
colnames(char_tab) = c("Social lineage", "Birds", "Songs" ,"Note classes")
xtable::xtable(char_tab, caption = "Table of the 10 social lineages, showing the number of individuals, the number of songs and the number of unique note classes found in all their songs.",
       label = "line_tab", digits = 0)

x$xtable