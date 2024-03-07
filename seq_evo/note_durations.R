#get duration data for optimal transport problem

head(unit_tab)

#we have 44.1kHz samples per second, windows of 512 length, hald overlapping

unit_tab = unit_tab %>%
  dplyr::mutate(nsam = floor(duration*44.1/256))

library(ggplot2)
ggplot(unit_tab, aes(x = note_label, y = nsam)) +
  geom_boxplot()

#do more exactly, using the note spectrograms directly
final_note_specgrams = unlist(final_note_specgrams, recursive = F)

dur_tab = lapply(final_note_specgrams, function(n){
  tibble::tibble(note_label = n$note_label, nsam = length(n$time))
})

dur_tab = do.call(rbind, dur_tab)

ggplot(dur_tab, aes(y = nsam, x = note_label)) +
  geom_boxplot()

#44kHz samples so...... 44*unit_tab$duration>512, i.e. duration cutoff is 512/44

# unit_tab = unit_tab %>%
#   dplyr::filter(duration > 512/44)