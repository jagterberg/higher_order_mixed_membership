###################################
#' This script runs all of the simulations and saves the outputs.  
#' It is used as part of producing figure 3 in the Appendix C.1 of the supplementary
#' material of the paper. It is run using tensor_two_infty.R, which runs one simulation 
#' for a given dimension p, noise level sigma_max, and amount of heteroskedasticity.
#' In this script we parallelize for different values of each of these and run the
#' simulations.
############################

source("misc.R")
source("HOOI.R")
source("membership_estimation")
source("tensor_two_infty.R")

##############################
# initialize parameters for simulations
numcores<- detectCores()
cl <- makeCluster(numcores,type = "FORK")
registerDoParallel(cl,numcores)

print("got here!")
rs <- rep(3,3)
set.seed(9192022)

#create core tensor
C <- array(rnorm(prod(rs)), dim=rs) 
C <- as.tensor(C)

# adjust minimal separation in core tensor
delta = 10
Sk = k_unfold(C, 1)@data
delta.min = min(svd(Sk)$d)
for (i in 1:3){
  Sk = k_unfold(C, i)@data
  for (k1 in 1:(rs[i]-1)){
    for (k2 in (k1+1):rs[i]){
      delta.min = min(delta.min, min(svd(Sk)$d))
    }
  }
}
C = C * delta / delta.min

#adjust ranks, dimensions, values of sigma, and number of simulations
r <- 3
ps <- seq(100,500,50)
ntrials <- 10
sigmas <- seq(1,100,5)

# we now run one simulation for each p in ps, for each number of trials, and for each
# different value of sigma above.

# for each p:
finalres_uniform <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  
  # for each number of trials:
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    
    # for each different noise value:
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      
      #run one simulation     
      return(runonesim(p,r,C,sigmamax = sigmas[j]))
     }
   return(toreturn)
  }
  return(toreturn2)
}
  
save(finalres_uniform,file = "../output/sim1_9-19.Rdata")
print("first sim done!")

# run again with more heteroskedasticity
finalres_het1 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .75))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het1,file = "../output/sim2_9-19.Rdata")
print("second sim done!")

# run again with even more heteroskedasticity:
finalres_het2 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .5))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het2,file = "../output/sim3_9-19.Rdata")
print("third sim done!")

# finally, include the most heteroskedasticity
finalres_het3 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .25))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het3,file = "../output/sim4_9-19.Rdata")
print("finished!")
