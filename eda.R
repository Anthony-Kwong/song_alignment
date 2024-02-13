#exploratory data analysis for sequences evo.

library(ggplot2)
head(unit_tab)
ggplot(data = unit_tab, aes(x = note_label, y = duration)) +
  geom_boxplot()
