pacman::p_load(stringr, collections, gplots)

setwd("./data/")

airquality
?filter()

months = unique(airquality$Month)
sapply(months, function(m){
  small_df = dplyr::filter(airquality, Month == m)
  mean(small_df$Wind)
})
