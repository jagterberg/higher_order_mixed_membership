#' This script has only one function: SPAMM, which stands for
#' "successive projection algorithm for mixed memberships."
#' It outputs three matrices of estimated mixed memberships.  
#' The first input is uhats, the list of estimated tensor singular vectors,
#' where each element of the list is a pk x rk orthonormal matrix.
#' The second input is a Boolean variable threshold, which provides the option
#' to threshold negative or extremely small values after solving for the estimated memberships.
#' the tval is the tolerance level.  In practice, having thresholding is useful
#' to ensure numerical stability.
SPAMM <- function(uhats,threshold=FALSE,tval=10^(-12)) {
  Pis <- list()
  Js <- list()
  
  # run the successive projection algorithm to find the vertices for each uhat
  for (k in c(1:length(uhats))) {
    R <- uhats[[k]]
    J <- NULL
    j = 1
    r <- dim(R)[2]
    while (j <= r) {
      jstar <- which.max(apply(R,1,function(x){
        sum(x^2)
      }))
      vj <- R[jstar,]
      R <- R %*% ( diag(1,r) - vj %*% t(vj) / sum(vj^2))
      J <- union(J,jstar)
      j <- j+1
    }
    
    #above finds the index set, we estimate the memberships
    # via the equation Pi = U %*% (U_J)^{-1}
    Pihat <- Uhats[[k]] %*% solve(Uhats[[k]][J,])
    
    # above often results in negative values and may not have unit ell_1 norms
    # so the next step is included to improve the estimation
    if(threshold) {
      #if the estimate is small or negative, we set the entry to zero
      Pihat <- ifelse(Pihat < tval,0,Pihat)
      
      # after this step, if all resulting entries are zero or very small, then
      # we may have issues in our normalization step. These estimates are not
      #very good regardless, so we will essentially ignore them for our normalization,
      # since we are attempting to divide zero by zero.  Effectively these are outliers,
      # and many other authors propose a pruning step to eliminate them.  We keep them
      # but we we will not output unit ell_1 norm estimates.
      
      # First we calculate row ell_1 norms:
      vec <- Pihat %*% rep(1,dim(Pihat)[2])
      
      # if the row norms are too small, then we set the row norms to be 1
      vec <- ifelse(vec < tval,1,vec)
      
      # we then normalize Pihat to have unit ell_1 norms
      Pihat <- diag(as.vector(1/vec),length(vec)) %*% Pihat
    } else {
      
      # if we don't threshold, then we are not concerned with negative or small
      # values:
      vec <- abs(Pihat) %*% rep(1,dim(Pihat)[2])
      Pihat <- diag(as.vector(1/vec),length(vec)) %*% Pihat
    }
    Pis[[k]] <- Pihat
    Js[[k]] <- J
  }
  return(list(Pis,Js))
}
