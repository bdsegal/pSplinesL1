
# compute confidence intervals for a fitted model
ci <- function(model, alpha = 0.05) {
  #' Obtain confidence intervals for l1ADMM objects
  #'
  #' This function obtains Bayesian and Frequentist confidence bands
  #' for models fit with the ADMMl1 function
  #' @param model output from l1ADMM with cv = FALSE
  #' @param alpha confidence level
  #' @keywords confidence intervals bands
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
  #' # run cross validation and view paths
  #' cvOut <- cv(y = "y", id = "id", X = X, rand = rand,
  #'                   K = 5,
  #'                   pathLength = 20,
  #'                   data = simData)
  #' cvOut
  #' 
  #' # fit model with all data
  #' a1 <- admm(y = "y", X, Z = rand$Z, S = rand$S,
  #'             tau = cvOut$smoothOpt[1],
  #'             lambda = cvOut$smoothOpt[2:(length(X)+1)],
  #'             rho = min(5, max(cvOut$smoothOpt)),
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

  J <- length(model$params$X)
  ord <- order(model$data$x)

  S <- model$params$S
  Z <- model$params$Z

  sigma2 <- model$fit$sigma2
  sigma2b <- sigma2 / model$params$tau

  V <- sigma2b * Z %*% tcrossprod(ginv(as.matrix(S)), Z) + diag(sigma2, nrow(Z))

  # smooths (beta_1 through beta_J), assume intercept should be added to beta_1 only
  smooth <- list()
  for (j in 1:J) {
    smooth[[j]] <- as.vector(model$param$X[[j]]$F %*% model$coef$beta[[j]])
  }
  smooth[[1]] <- smooth[[1]] + model$coefs$beta0
  
  keep <- list()
  yLowerBayesQuick <- list()
  yUpperBayesQuick <- list()
  yUpperPoint <- list()
  yLowerPoint <- list()

  for (j in 1:J) {
    print(j)
    
    F <- model$params$X[[j]]$F
    D <- model$params$X[[j]]$D
    Lambda <- model$params$lambda[j] * crossprod(D, D)

    if(j == 1) {
      F <- cbind(1, F)
      D <- bdiag(0, D)
      Lambda <- bdiag(0, Lambda)
    }

    keep[[j]] <- which(rowSums(F) != 0)

    W <- crossprod(F, solve(V, F)) + Lambda
    Winv <- ginv(as.matrix(W)) # ginv probably not necessary, but safer coding
    HtHnormDiagBayes <- sqrt(rowSums((F %*% Winv) * F))
   
    H <- F %*% solve(crossprod(F, F) + Lambda, t(F))
    HtHnormDiag <- sqrt(rowSums((H %*% V) * H))

    zQuant <- qnorm(1-alpha/2)
    yUpperBayesQuick[[j]] <- as.vector(smooth[[j]] + zQuant * HtHnormDiagBayes)
    yLowerBayesQuick[[j]] <- as.vector(smooth[[j]] - zQuant * HtHnormDiagBayes)
    yUpperPoint[[j]] <- as.vector(smooth[[j]] + zQuant * HtHnormDiag)
    yLowerPoint[[j]] <- as.vector(smooth[[j]] - zQuant * HtHnormDiag)
  }

  CI <- list()
  for (j in 1:J) {
    ordKeep <- ord[which(ord %in% keep[[j]])]
    CI[[j]] <- data.frame(x =  model$data$x[ordKeep],
                          smooth = smooth[[j]][ordKeep],
                          yUpperBayesQuick = yUpperBayesQuick[[j]][ordKeep],
                          yLowerBayesQuick = yLowerBayesQuick[[j]][ordKeep],
                          yLowerPoint = yLowerPoint[[j]][ordKeep],
                          yUpperPoint = yUpperPoint[[j]][ordKeep]
                          )
  }

  class(CI) <- "ci"

  return(CI)
}

print.ci <- function(CI) {
  #' Print function for confidence bands
  #'
  #' This function plots the confidence bands
  #' @param CI output from the ci function
  #' @keywords ci print
  #' @export
  J <- length(CI)
  for (j in 1:J) {
    CIpoly <- data.frame(x = c(CI[[j]]$x, rev(CI[[j]]$x)), 
                         y = c(CI[[j]]$yLowerBayesQuick, 
                               rev(CI[[j]]$yUpperBayesQuick)))

    yRange <- range(c(CI[[j]]$yLowerBayesQuick, 
                       CI[[j]]$yUpperBayesQuick))

    par(mar = c(5, 5, 4, 2))
    plot(x = CI[[j]]$x, y = CI[[j]]$smooth, type = "n",
         ylim = yRange,
         xlab = "x", ylab = "y", main = paste("Smooth", j),
         cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.6)
    polygon(CIpoly$x, y = CIpoly$y, col = "gray", border = NA)
    lines(x = CI[[j]]$x, y = CI[[j]]$smooth, lwd = 2)

    readline("Hit any key to see next plot")
  }
}