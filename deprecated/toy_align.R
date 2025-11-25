#converting Mark's alignment code to R and adding comments

#we need to build a set of functions

#logsum function----

#use the log-sum-exp trick to compute log likelihoods without overflow/underflow

#input: numeric vector of likelihoods
#output: sum of the loglikelihoods in input
logsum <- function(v){
  vmax = max(v)
  res = vmax + log(sum( exp(v - vmax) ))
  return(res)
}

#test 
v = c(2000,2000)
log(sum(v)) #overflow
logsum(v) #correct

v = c(-2000,-2000)
log(sum(v)) #underflow
logsum(v) #correct

#other stuff
#not the sort of alignment tool I was looking for but could be useful later
library(DescTools)
x = c("AABBC", "ABBBC", "BBC")
cbind(StrAlign(x, sep = "\\l"))


