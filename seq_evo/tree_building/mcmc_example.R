#script with worked example for getting likelihood of observing one aligned (father + son) position, 
#using the Hamiltonian MCMC method

#there are two note classes, A and B

#father's counts at the position, 5A's and 5B's
N = c(5,5)
#son's counts, 7A's and 3B's
M = c(7,3)

#transition matrix T, (made up for now), we have estimated a proper one on the cluster
T = matrix(c(0.9,0.1,0.1,0.9), nrow = 2)



