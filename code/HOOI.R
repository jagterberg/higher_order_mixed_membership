###############################################
#' This script has two functions: HOOI_dd and HOOI_vanilla.  HOOI_dd runs HOOI
#' with the diagonal-deletion initialization as proposed in the main paper.
#' HOOI_vanilla runs HOOI with either a user-provided initialization or the SVD
#' initialization (which is the default).  HOOI_dd just computes the initialization
#' and then passes these initializations to HOOI_vanilla.
#################################################


#' Function to run diagonal-delted HOOI.  
#' Inputs: tens -- an order three tensor, 
#' r -- a vector of three integers of ranks
#' niter -- maximum number of iterations to run
#' output: estimates of the tensor singular vectors
HOOI_dd <- function(tens,r=c(1,1,1),niter=NULL) {
  # first, calculate the gram matrix of first matricization and set diagonal to zero
  to_svd <- k_unfold(tens,1)@data
  to_svd <- to_svd %*% t(to_svd)
  diag(to_svd) <- 0
  
  # next, compute the leading r[1] left singular vectors of this obtained matrix
  init1 <- irlba(to_svd,r[1])$u
  
  # repeat for modes two and three:
  to_svd <- k_unfold(tens,2)@data
  to_svd <- to_svd %*% t(to_svd)
  diag(to_svd) <- 0
  init2 <- irlba(to_svd,r[2])$u
  
  #with these initializations, pass them two HOOI_vanilla:
  return(HOOI_vanilla(tens=tens,r=r,niter=niter,init1=init1,init2=init2))
}


#' Function to run HOOI with either use-provided initialization or the SVD
#' initialization (default if init1, init2 are both null).  
HOOI_vanilla <- function(tens,r=c(1,1,1),niter=NULL,init1=NULL,init2=NULL) {
  
  #obtain the dimensions of each mode
  p1 <- dim(tens)[1]
  p2 <- dim(tens)[2]
  p3 <- dim(tens)[3]
  r1 <- r[1]
  r2 <- r[2]
  r3 <- r[3]
  
  #set the number of iterations if not specified
  if(is.null(niter)) {
    niter <- max(5,log(max(p1,p2,p3))) #uses logarithmically many iterations
  }
  
  #get initializations if not provided
  if(is.null(init1)) {
    print("obtaining initialization 1 by vanilla (truncated) SVD...")
    init1 <- irlba(k_unfold(tens,1)@data,r1)$u
  } 
  
  if(is.null(init2)) {
    print("obtaining initialization 1 by vanilla (truncated) SVD...")
    init2 <- irlba(k_unfold(tens,2)@data,r2)$u
  }
  
  print("beginning power iteration...")
  
  #using these two initializations, obtain Uhat3, the initialization for the third
  # mode
  to_irlba <- ttm( ttm(tens,t(init1),m=1),t(init2),m=2)
  dim3 <- k_unfold(to_irlba,3)@data
  Uhat3 <- svd(dim3,r3)$u
  Uhat1 <- init1
  Uhat2 <- init2
  i <- 1
  
  #now, we iterate by running HOOI for niter iterations
  while (i < niter) {
    print(paste("iteration",i,"of",niter))
    to_irlba <- ttm(ttm(tens,t(Uhat3),m=3),t(Uhat2),m=2)
    dim1 <- k_unfold(to_irlba,1)@data
    Uhat1 <- svd(dim1)$u[,c(1:r1)]
    
    to_irlba <- ttm(ttm(tens,t(Uhat1),m=1),t(Uhat3),m=3)
    dim2 <- k_unfold(to_irlba,2)@data
    Uhat2 <- svd(dim2)$u[,c(1:r2)]
    
    to_irlba <- ttm(ttm(tens,t(Uhat1),m=1),t(Uhat2),m=2)
    dim3 <- k_unfold(to_irlba,3)@data
    Uhat3 <- svd(dim3)$u[,c(1:r3)]
    i <- i+1
  }
  
  # return a list of estimated singular vectors
  return(list(Uhat1,Uhat2,Uhat3))
}
