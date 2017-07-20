
soft <- function(kappa, x) {
  #' soft thresholding function
  #'
  #' @param kappa scalar threshold
  #' @param x scalar value to be soft thresholded
  #' @keywords soft thresholding
  #' @export
  
  (x - kappa) * (x > kappa) +
  0 * (abs(x) <= kappa) +
  (x + kappa) * (x < -kappa)
}

msub <- function(i, p) {
  #' sub function for making inverse of finite difference matrix
  #'
  #' @keywords Finite difference inverse
  #' @export

  Ltemp <- matrix(1, nrow = p - i, ncol = p - i)
  Ltemp[upper.tri(Ltemp)] <- 0
  bdiag(diag(i), Ltemp)
}

Mfun <- function(k, p) {
  #' Function for making inverse of finite difference matrix
  #'
  #' @param k scalar k+1 is order of finite difference matrix
  #' @param p scalar column dimension
  #' @keywords Finite difference inverse
  #' @export
  
  mTemp <- msub(0, p) 
  for (i in 1:k) {
    mTemp <- mTemp %*%  msub(i, p)
  }
  return(mTemp)
}
