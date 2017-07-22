
# with random effects,  -------------------------------------------------------
admm <- function(y, 
                 X, 
                 Z,
                 S,
                 lambda,
                 tau,
                 data, 
                 rho = 1,
                 epsilonAbs = 1e-4,
                 epsilonRel = 1e-4,
                 iterMax = 1e3,
                 warm = NULL,
                 forCV = FALSE,
                 centerZ = FALSE, 
                 verbose = FALSE) {
  #' ADMM algorithm for fitting l1-penalized additive mixed models
  #'
  #' This function obtains point estimates via the ADMM algorithm
  #' @param y character string denoting output variable
  #' @param X list of fixed effect design matrices setup by l1smooth
  #' @param Z random effect design matrix setup by randomEffects
  #' @param S random effect penalty matrix setup by randomEffects
  #' @param lambda scalar smoothing parameter for random effects
  #' @param tau vector of smoothing parameters for fixed effects
  #' @param rho ADMM penalty parameter
  #' @param epsilonAbs absolute error in ADMM (tolerance for convergence)
  #' @param epsilonRel relative error in ADMM (tolerance for convergence)
  #' @param iterMAx maximum number of ADMM iterations
  #' @param warm list of warm start parameters
  #' @param data dataset
  #' @param forCV TRUE/FALSE If FALSE, returns less information in fitted object
  #' @param centerZ TRUE/FALSE If TRUE, imposes centering constraints on random effects
  #' @param verbose TRUE/FALSE If TRUE, prints ADMM iteration
  #' @keywords ADMM
  #' @export
  #' @examples
  #' library(psplinesl1)
  #' data(simData)
  #'
  #' # setup p-spline matrices
  #' X <- list(ps(x = "x", data = simData, 
  #'              norder = 2, k = 1, width = 0.05,
  #'              center = TRUE))
  #'
  #' # setup random effect matrices
  #' rand <- re(x = "x", id = "id", data = simData,
  #'            randomCurves = FALSE)
  #'
  #' # fit model with ADMM
  #' a1 <- admm(y = "y", X, Z = rand$Z, S = rand$S,
  #'             tau = 0.0502,
  #'             lambda = 0.169,
  #'             rho = 0.169,
  #'             data = simData)
  #'
  #' # get and plot fitted model with confidence bands
  #' CI <- ci(model = a1, alpha = 0.05)
  #' CI
  #' 
  #' # extract values from ci object for custom plotting
  #' CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
  #'                      y = c(CI[[1]]$yLowerBayesQuick, 
  #'                            rev(CI[[1]]$yUpperBayesQuick)))
  #' 
  #' ggplot(aes(x = x, y = y), data = newDat)+
  #'   geom_polygon(data = CIpoly, fill = "grey")+
  #'   geom_line(aes(x = CI[[1]]$x, y = CI[[1]]$smooth))

  ycol <- which(colnames(data) == y)

  if(length(X) != length(lambda)) {
    warning("number of smooths must equal number of smoothing parameters (lambda)")
  }

  J <- length(X)

  # centering constraints
  if(centerZ) {
    Q <- qr.Q(qr(colSums(Z)), complete = TRUE)[, -1]
    Z <- Z %*% Q
    S <- crossprod(Q, S %*% Q)
  }

  Finter <- list()
  Fqr <- list()
  for (j in 1:J) {
    Finter[[j]] <- crossprod(X[[j]]$F) + rho * crossprod(X[[j]]$D)
    qrDecomp <- qr(Finter[[j]], tol = 1e-10)
    Fqr[[j]] <- list(Qt = t(qr.Q(qrDecomp)),
                           R = qr.R(qrDecomp))
  }

  qrDecomp <- qr(crossprod(Z) + tau * S, tol = 1e-10)
  if ("sparseQR" %in% class(qrDecomp)) {
    ZSqr <- list(Qt = t(qr.Q(qrDecomp)), 
                 R = qrR(qrDecomp))
  } else {
    ZSqr <- list(Qt = t(qr.Q(qrDecomp)), 
                 R = qr.R(qrDecomp))
  }

  if (is.null(warm)) {
    u <- list()
    wNew <- list()
    for (j in 1:J) {
      u[[j]] <- rep(0, X[[j]]$p - X[[j]]$k - 1)
      wNew[[j]] <- rep(0, X[[j]]$p - X[[j]]$k - 1)
    }
    beta <- list()
    for (j in 1:J) {
      beta[[j]] <- rep(0, ncol(X[[j]]$F))
    }
  } else {
    u <- warm$u
    wNew <- warm$wNew
    beta <- warm$beta
  }
  b <- rep(0, ncol(Z))

  beta0 <- 0
  w <- list()

  conv <- FALSE
  iter <- 1
  rNorm <- rep(NA, iterMax)
  sNorm <- rep(NA, iterMax)

  dimPrim <- sum(sapply(X, function(x) {x$p - x$k})) - J
  dimDual <- sum(sapply(X, function(x) {x$p}))

  while(!conv & iter <= iterMax) {
    if(verbose){print(iter)}

    # beta0 update
    yRes <- data[, ycol] - Z %*% b
    for (j in 1:J) {
      yRes <- yRes - X[[j]]$F %*% beta[[j]]
    }
    beta0 <- mean(yRes)
    
    # move wNew to w
    for(j in 1:J) {
      w[[j]] <- wNew[[j]]
    }

    # beta updates
    for (j in 1:J) {
      yRes <- data[, ycol] - beta0 - Z%*%b
        for (l in 1:J) {
          if (l != j) {
            yRes <- yRes - X[[l]]$F %*% beta[[l]]
          }
        }
      beta[[j]] <- backsolve(Fqr[[j]]$R, Fqr[[j]]$Qt %*% (crossprod(X[[j]]$F, yRes) + 
                             rho * crossprod(X[[j]]$D, w[[j]] - u[[j]])))
    }

    # w updates
    for (j in 1:J) {
      wNew[[j]] <- soft(lambda[j] / rho, X[[j]]$D %*% beta[[j]] + u[[j]])
    }

    # primal residuals
    rList <- list()
    for (j in 1:J) {
      rList[[j]] <- as.vector(X[[j]]$D %*% beta[[j]] - wNew[[j]])
    }
  
    # u update
    for (j in 1:J) {
      u[[j]] <- as.vector(u[[j]] + rList[[j]])
    }

    # dual residuals
    sList <- list()
    for (j in 1:J) {
      sList[[j]] <- as.vector(-rho * crossprod(X[[j]]$D, wNew[[j]] - w[[j]]))
    }

    # random effects update
    yRes <- data[, ycol] - beta0
    for (j in 1:J) {
      yRes <- yRes - X[[j]]$F %*% beta[[j]]
    }

    # update b
    b <- backsolve(ZSqr$R, ZSqr$Qt %*% crossprod(Z, yRes))

    # get primal and dual residuals
    r <- do.call(c, rList)
    s <- do.call(c, sList)

    Dbeta <- list()
    for (j in 1:J) {
      Dbeta[[j]] <- X[[j]]$D %*% beta[[j]]
    }
    DbetaBind <- do.call(rbind, Dbeta)
    wBind <- do.call(rbind, wNew)

    epsilonPri <- sqrt(dimPrim) * epsilonAbs + 
                  epsilonRel * max(norm(DbetaBind, type = "F"),
                                   norm(wBind, type = "F"))

    Du <- list()
    for (j in 1:J) {
      Du[[j]] <- crossprod(X[[j]]$D, u[[j]])
    }
    DuBind <- do.call(rbind, Du)
    epsilonDual <- sqrt(dimDual) * epsilonAbs + 
                   epsilonRel * rho * norm(DuBind, type = "F")                   

    # iterate until true
    rNorm[iter] <- as.numeric(sqrt(crossprod(r)))
    sNorm[iter] <- as.numeric(sqrt(crossprod(s)))
    conv <- (rNorm[iter] < epsilonPri) & (sNorm[iter] < epsilonDual)
    iter <- iter + 1
  }

  yMarg <- beta0
  for (j in 1:J) {
    yMarg <- yMarg + X[[j]]$F %*% beta[[j]]
  }
  yMarg <- as.vector(yMarg)
  yHat <- as.vector(yMarg + Z %*% b)
  residuals <- as.vector(data[, ycol] - yHat)
  
  if(forCV) {
    return(list(conv = list(converged = conv,
                            iter = iter - 1,
                            rNorm = rNorm[1:(iter - 1)],
                            sNorm = sNorm[1:(iter - 1)],
                            epsilonAbs = epsilonAbs,
                            epsilonRel = epsilonRel),
                warm = list(u = u,
                            wNew = wNew,
                            beta = beta,
                            beta0=beta0),
                fit = list(residuals = residuals)
              )
          )
  } else {

    SSE <- sum(residuals^2)
    n <- table(data$id)
    sumn <- sum(n)

    dfTable <- dfl1(tau = tau, lambda = lambda, X = X, wNew = wNew, Z = Z, S = S)
   
    # use Stein's method if possible, and restricted approximation if necessary
    if(!is.na(dfTable$Overall[which(rownames(dfTable) == "Stein")])) {
      dfSSE <- sumn - dfTable$Overall[which(rownames(dfTable) == "Stein")]
    } else {
      dfSSE <- sumn - dfTable$Overall[which(rownames(dfTable) == "Restricted")]
    }

    # use ridge approximation if possible, and restricted approximation if necessary
    if(!is.na(dfTable$Overall[which(rownames(dfTable) == "Ridge")])) {
      dfSSEridge <- sumn - dfTable$Overall[which(rownames(dfTable) == "Ridge")]
    } else {
      dfSSEridge <- sumn - dfTable$Overall[which(rownames(dfTable) == "Ridge restricted")]
    }

    sigma2 <- SSE / dfSSE
    sigma2Ridge <- SSE / dfSSEridge

    return(list(conv = list(converged = conv,
                            iter = iter - 1,
                            rNorm = rNorm[1:(iter - 1)],
                            sNorm = sNorm[1:(iter - 1)],
                            epsilonAbs = epsilonAbs,
                            epsilonRel = epsilonRel),
                data = data,
                coefs = list(beta0 = beta0,
                             beta = beta,
                             b = b),
                params = list(lambda = lambda,
                              rho = rho,
                              tau = tau,
                              X = X,
                              Z = Z,
                              S = S,
                              u = u,
                              wNew = wNew),
                fit = list(residuals = residuals,
                           yHat = yHat,
                           df = dfTable,
                           sigma2 = sigma2,
                           sigma2Ridge = sigma2Ridge,
                           dfSSE = dfSSE))
          )
  }
}
