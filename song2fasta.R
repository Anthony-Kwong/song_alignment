#convert song sequences to fasta (example)
b37_a

alignment_37
a_fasta = bio3d::as.fasta(alignment_37)
bio3d::write.fasta(a_fasta, file = "./data/songs_fasta/JS0037.fasta")
ggmsa("./data/songs_fasta/JS0037.fasta", color = "LETTER") + geom_seqlogo() + geom_msaBar()

library(seqinr)
write.fasta(b37_a, names = paste("JS0037", seq(1:10),sep="_"), file.out = "./data/songs_fasta/JS0037.fasta", as.string = F)

song_sequences <- system.file("./data/songs_fasta/JS0037.fasta")
ggmsa("./data/songs_fasta/JS0037.fasta") + geom_seqlogo() + geom_msaBar()

protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
ggmsa(protein_sequences, start = 221, end = 280, char_width = 0.5, seq_name = T) + geom_seqlogo() + geom_msaBar()