#check stats

table(bird_songs$Line)

#table of song length statistics for every lineage ----
#get song length statistics for every lineage

song_stats = bird_songs %>%
  dplyr::group_by(Line) %>%
  dplyr::summarise(
    `median length` = median(nchar(note.seq), na.rm = TRUE),
    `minimum length`    = min(nchar(note.seq), na.rm = TRUE),
    `maximum length`    = max(nchar(note.seq), na.rm = TRUE),
    .groups = "drop"
  )

xtable::xtable(song_stats, caption = "Table of song length statistics for every lineage.",
               label = "tab:songlen_tab", digits = 0)
