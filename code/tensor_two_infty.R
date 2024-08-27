############################################
#' This script contains a single function that runs one single
#' simulation for the simulations in Appendix C.1 in the supplementary material.
#' This is part of the code used to generate Figure 3.
############################################

#' input: dimensions p, ranks r, core tensor C, maximum standard deviation sigmamax,
#' and amount of heteroskedasticity het_amt.
#' Output: ell_{2,infty} error between true and estimated memberships (minimized over
#' permutations)
runonesim <- function(p,r,C,sigmamax,het_amt = 1) {
  #generate pure nodes
  purenodes <- c(1,rep(0,r-1))
  for(i in c(2:r)) {
    newnode <- c( rep(0,(i-1)),1,rep(0,(r-i)))
    purenodes <- rbind(purenodes,newnode)
  }
  
  #generate memberships from a dirichlet distribution
  Pi1 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  Pi2 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  Pi3 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  
  #true tensor: underlying r x r x r rank tensor
  T_true <- ttm(ttm(ttm(C,Pi1,m=1),Pi2,m=2),Pi3,m=3)
  
  # generate standard deviations uniformly if no significant heteroskedasticity
  if (het_amt == 1) {
    sds <- runif(p^3,0,sigmamax)
  } else {
    # otherwise generate them from a beta distribution
    sds <- sigmamax*rbeta(p^3,het_amt,het_amt)
  }
  
  # generate noise tensor
  Z <- as.tensor(array(rnorm(p^3,sd=sds),dim=c(p,p,p)))
  T_obs <- T_true + Z
 
  # estimate memberships 
  clusts <- SPAMM(uhats = HOOI_dd(T_obs,r=c(r,r,r)))
  
  #take absolute values
  Pihat1 <- abs(clusts[[1]][[1]])
  
  #minimize error over permutations
  perms <- rbind(
    c(1,2,3),
    c(1,3,2),
    c(2,1,3),
    c(2,3,1),
    c(3,2,1),
    c(3,1,2)
  )
  minvals <- rep(0,nrow(perms))
  for ( i in c(1:nrow(perms))) {
    minvals[i] <- sqrt(sum( (Pi1-Pihat1[,perms[i,]])^2 ))
  }
  
  return(min(minvals))
  
}



