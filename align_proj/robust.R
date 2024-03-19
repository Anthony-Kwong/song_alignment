#investigate robustness of alignments to different parameters

bird_songs = readr::read_csv("~/Documents/GitHub/song_alignment/data/NoteSequences.csv")
birds = bird_songs$Bird.ID

js = "JS0329"

birds = c("JS0329","JS0037")

js_songs = bird_songs %>%
  dplyr::filter(Bird.ID == js)

setwd("~/Documents/GitHub/song_alignment/results/fasta")

purrr::walk(.x, .f)

#testing parameters----


#map method for gaps
lambda = c(0,5,10)

#go through every bird, generate PHMM alignment with different lambdas and plot
purrr::walk(birds, function(js){
  #set model parameters to test----
  
  #filter to bird of interest's songs
  js_songs = bird_songs %>%
    dplyr::filter(Bird.ID == js)
  
  #preprocessing----
  
  #get note sequences as long strings
  bird_songsseqs = js_songs$note.seq
  
  #get sequences in split strings
  bird_songs_split = lapply(bird_songsseqs, function(s){strsplit(s, "")[[1]]})
  bird_songs_split <- bird_songs_split[order(lengths(bird_songs_split),decreasing = T)]
  
  #get the alphabet
  comb = paste(bird_songsseqs, collapse = "")
  letters = unique(strsplit(comb, "")[[1]])
  
  ###alignment with different parameters
  
  # tidyr::crossing(lambda, opt_method)
  
  #loop through the 2 optimization procedures
  
  #fit the models for different lambdas
  bw_models = lapply(lambda, function(l){
    derivePHMM(bird_songs_split, residues = letters, pseudocounts = "Laplace", refine = "BaumWelch",
               inserts = "map", lambda = l)
  })
  
  #create a new directory for the bird if it doesn't exist already
  bird_dir = paste("./robust_test/map/",js, sep = "")
  dir.create(bird_dir, showWarnings = FALSE)
  
  #loop thru lambdas and generate plots
  for(i in 1:length(lambda)){
    alignment = align(bird_songs_split, model = map_models[[i]], seqweights = NULL, residues = letters)
    A = og_order(align_mat = alignment, song_seqs = bird_songsseqs)
    alignment_fasta = bio3d::as.fasta(A)
    fname = paste(bird_dir,"/", js,".fasta",sep="")
    print(fname)
    bio3d::write.fasta(alignment_fasta, file = fname)
    
    #make plot
    plot = ggmsa(fname, color = "LETTER") + geom_msaBar() 
    plotname = paste(bird_dir,"/",js,"_",lambda[i],".png", sep = "")
    ggsave(plot, file = plotname)
  }
})
