#######################
#' This script simply loads the following packages: irlba,
#' Matrix, rTensor, gtools, foreach, and doParallel.  It also loads two functions: 
#' procrustes and two_infty. procrustes calculates the orthogonal matrix 
#' optimally aligning two orthonormal matrices with respect to the Frobenius 
#' norm, and two_infty calculates the ell_two_infty distance
#' between two orthonormal matricecs after Procrustes rotation.
###########################

#Load packages
if(!require(irlba)) {
  install.packages("irlba")
  library(irlba)
}
if(!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}
if(!require(rTensor)) {
  install.packages("rTensor")
  library(rTensor)
}
if(!require(gtools)) {
  install.packages("gtools")
  library(gtools)
}
if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(doParallel)) {
  install.packages("doParallel")
  library(doParallel)
}

# calculate the optimal rotation W = argmin_{W} || Uhat W - U ||_F, which
# can be computed from the svd of Uhat^T U as seen below
procrustes <- function(Uhat,U) {
  inner_prod <- t(Uhat) %*% U
  toReturn <- svd(inner_prod)
  
  return(toReturn$u %*% toReturn$v)
}

# This function calculates the (approximate) ell two infty distance between
# two orthonormal matrices.  In lieu of the optimal ell two infty orthogonal matrix
# it uses the optimal Frobenius orthogonal matrix from procrustes().  
two_infty <- function(Uhat,U) {
  W <- procrustes(Uhat,U)
  norms <- apply(Uhat %*% W-U,1,function(x){
    sum(x^2)
  })
  return(
    max(sqrt(norms))
  )
}





