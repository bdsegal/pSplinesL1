dfl1 <- function(tau, lambda, X, wNew, Z, S) {
  #' Obtain degrees of freedom for l1ADMM objects
  #'
  #' This function obtains degrees of freedom estimates
  #' for models fit with the admm function. Usually called from within admm.
  #' @param tau vector of smoothing parameters for fixed effects
  #' @param lambda scalar smoothing parameter for random effects
  #' @param X list of fixed effect design matrices setup by l1smooth
  #' @param wNew list of updated w parameters in ADMM algorithm
  #' @param Z random effect design matrix setup by randomEffects
  #' @param S random effect penalty matrix setup by randomEffects
  #' @keywords confidence intervals bands
  #' @export
  #' @examples
  #' dfTable <- dfl1(tau, lambda, X, wNew, Z, S)

  J <- length(X)

  # degrees of freedom --------------------------------------------------------
  activeSet <- list()
  for (j in 1:J) {
    activeSet[[j]] <- c(1:(X[[j]]$k+1), (X[[j]]$k+1) + which(wNew[[j]] != 0))
  }
  asLen <- sapply(activeSet, length)
  
  H <- list()
  Vj <- list()
  for (j in 1:J) {
    H[[j]] <- X[[j]]$F %*% tcrossprod(ginv(as.matrix(crossprod(X[[j]]$F) + 
                         lambda[j] * crossprod(X[[j]]$D))), X[[j]]$F)
    M <- Mfun(k = X[[j]]$k, p = ncol(X[[j]]$F))
    Vj[[j]] <- as.matrix((X[[j]]$F %*% M)[, activeSet[[j]]])
  }

  Hb <- Z %*% solve(crossprod(Z) + tau * S, t(Z), sparse = TRUE)

  A <- cbind(do.call(cbind, Vj), Z)
  activeSetCum <- c(0, cumsum(asLen), ncol(A))

  # Stein degrees of freedom -------------------------------
  Omega <- array(0, dim = rep(ncol(A), 2))
  Omega[((activeSetCum[J+1]+1):activeSetCum[J+2]), 
        ((activeSetCum[J+1]+1):activeSetCum[J+2])]  <- tau * as.matrix(S)

  dfStein <- NA
  dfSteinj = rep(NA, J+1)
  try({
    AtA <- crossprod(A)
    Adiag <- diag(solve(AtA + Omega, AtA))
    for (j in 1:(J+1)) {
      dfSteinj[j] <- sum(Adiag[(activeSetCum[j]+1):activeSetCum[j+1]])
    }
    dfStein <- sum(dfSteinj) + 1
  })

  # ADMM approximation ------------------------------------
  kVec <- sapply(1:length(X), function(j) {X[[j]]$k})
  byVec <- as.numeric(sapply(1:length(X), function(j) {is.null(X[[j]]$Q)}))
  dfADMMj <- c(sapply(wNew, function(x) sum(x != 0)) + kVec + byVec,
              sum(diag(Hb)))
  dfADMM <- sum(dfADMMj) + 1

  # restricted derivatives approximation ------------------
  dfResj <- c(sapply(Vj, function(V) {sum(solve(crossprod(V), crossprod(V)))}), sum(diag(Hb)))
  dfRes <- 1 + sum(dfResj)

  # ridge approximation ----------------------------------
  Flist <- list()
  for (j in 1:J) {
    Flist[[j]] <- X[[j]]$F
  }
  U <- cbind(do.call(cbind, Flist), Z)

  p <- c(sapply(Flist, ncol), ncol(Z))
  pCum <- c(0, cumsum(p))

  OmegaRidge <- array(0, dim = rep(ncol(U), 2))
  for (j in 1:J) {
    OmegaRidge[((pCum[j]+1):pCum[j+1]), ((pCum[j]+1):pCum[j+1])] <- lambda[j] * 
      crossprod(X[[j]]$D)
  }
  OmegaRidge[((pCum[J+1]+1):pCum[J+2]),((pCum[J+1]+1):pCum[J+2])]  <- tau * as.matrix(S)

  # true ridge approximation
  dfRidge <- NA
  dfRidgej = rep(NA, J+1)
  try({
    UtU <- crossprod(U)
    Hdiag <- diag(solve(UtU + OmegaRidge, UtU))
    for (j in 1:(J+1)) {
      dfRidgej[j] <- sum(Hdiag[(pCum[j]+1):pCum[j+1]])
    }
    names(dfRidgej) <- c(paste("F", 1:length(X), sep = ""), "Z")
    dfRidge <- sum(dfRidgej) + 1
  })

  # restricted derivative ridge approximation -------------
  dfRidgeResj <- c(sapply(H, function(h) {sum(diag(h))}), sum(diag(Hb)))
  dfRidgeRes <- 1 + sum(dfRidgeResj)
  # dfRidgeResLambda <- sum(dfRidgeResj[-length(dfRidgeResj)])

  dfTable <- data.frame(rbind(c(dfStein, dfSteinj),
                              c(dfRes, dfResj),
                              c(dfADMM, dfADMMj),
                              c(dfRidge, dfRidgej),
                              c(dfRidgeRes, dfRidgeResj)))
  colnames(dfTable) <- c("Overall", c(paste("F", 1:J, sep = ""), "Z"))
  rownames(dfTable) <- c("Stein", "Restricted", "ADMM", "Ridge", "Ridge restricted")

  return(dfTable)
}
