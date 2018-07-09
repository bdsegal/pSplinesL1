ci <- function(model, alpha = 0.05, newData = NULL, lenOut = 100) {
  #' Obtain confidence intervals for l1ADMM objects
  #'
  #' This function obtains Bayesian and Frequentist confidence bands
  #' for models fit with the ADMMl1 function
  #' @param model output from l1ADMM with cv = FALSE
  #' @param alpha confidence level
  #' @param newData new dataset at which to get the predicted values and confidence intervals
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
  #' plot(cvOut)
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
  #' plot(CI)
  #' 
  #' # extract values from ci object for custom plotting
  #' CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
  #'                      y = c(CI[[1]]$lower, 
  #'                            rev(CI[[1]]$upper)))
  #' 
  #' ggplot(aes(x = x, y = y), data = newDat)+
  #'   geom_polygon(data = CIpoly, fill = "grey")+
  #'   geom_line(aes(x = CI[[1]]$x, y = CI[[1]]$smooth))

  covariates <- unique(do.call(c, lapply(model$params$X, function(x){c(x$x, x$by)})))

  if (!is.null(newData) & any(!covariates %in% colnames(newData))) {
    stop("New dataset must contain all predictors in original dataset")
  }

  J <- length(model$params$X)
  Xnew <- list()
  xSorted <- list()

  for (j in 1:J) {
    if (!is.null(newData)) {
      xcol <- which(colnames(newData) == model$params$X[[j]]$x)
      xSorted[[j]] <- sort(unique(newData[, xcol]))
    } else {
      xcol <- which(colnames(model$data) == model$params$X[[j]]$x)
      xSorted[[j]] <- sort(unique(model$data[, xcol]))
    }

    predictData <- data.frame(x = xSorted[[j]])

    # always show the effect of a 1 unit increase in by variable
    byVal <- model$params$X[[j]]$by
    predictData[, byVal] <- 1

    Xnew[[j]] <- psSub(xcol, 
                       x = model$params$X[[j]]$x,
                       basis = model$params$X[[j]]$basis, 
                       norder = model$params$X[[j]]$norder,
                       k = model$params$X[[j]]$k, 
                       data = predictData,
                       by = model$params$X[[j]]$by,
                       ref = model$params$X[[j]]$ref,
                       width = model$params$X[[j]]$width,
                       center = model$params$X[[j]]$center)
  }

  sigma2 <- model$fit$sigma2
  if (is.null(model$fit$sigma2b)) {
    sigma2b <- sigma2 / model$params$tau
  } else {
    sigma2b <- model$fit$sigma2b
  }
  
  S <- model$params$rand$S
  Z <- model$params$rand$Z 
  V <- sigma2b * Z %*% tcrossprod(ginv(as.matrix(S)), Z) + diag(sigma2, nrow(Z))

  # smooths (beta_1 through beta_J), assume intercept should be added to beta_1 only
  smooth <- list()
  for (j in 1:J) {
    smooth[[j]] <- as.vector(Xnew[[j]]$F %*% model$coef$beta[[j]])
  }
  smooth[[1]] <- smooth[[1]] + model$coefs$beta0
  
  yLowerBayesQuick <- list()
  yUpperBayesQuick <- list()
  yUpperPoint <- list()
  yLowerPoint <- list()

  for (j in 1:J) {

    Fnew <- Xnew[[j]]$F
    F <- model$params$X[[j]]$F
    D <- model$params$X[[j]]$D
    Lambda <- model$params$lambda[j] * crossprod(D)

    if(j == 1) {
      Fnew <- cbind(1, Fnew)
      F <- cbind(1, F)
      D <- bdiag(0, D)
      Lambda <- bdiag(0, Lambda)
    }

    W <- crossprod(F, solve(V, F)) + Lambda
    Winv <- ginv(as.matrix(W)) # ginv probably not necessary, but safer coding
    HtHnormDiagBayes <- sqrt(rowSums((Fnew %*% Winv) * Fnew))
   
    H <- Fnew %*% solve(crossprod(F) + Lambda, t(F))
    HtHnormDiag <- sqrt(rowSums((H %*% V) * H))

    zQuant <- qnorm(1-alpha/2)
    yUpperBayesQuick[[j]] <- as.vector(smooth[[j]] + zQuant * HtHnormDiagBayes)
    yLowerBayesQuick[[j]] <- as.vector(smooth[[j]] - zQuant * HtHnormDiagBayes)
    yUpperPoint[[j]] <- as.vector(smooth[[j]] + zQuant * HtHnormDiag)
    yLowerPoint[[j]] <- as.vector(smooth[[j]] - zQuant * HtHnormDiag)
  }

  CI <- list()
  for (j in 1:J) {
    CI[[j]] <- data.frame(x =  xSorted[[j]],
                          smooth = smooth[[j]],
                          lower = yLowerBayesQuick[[j]],
                          upper = yUpperBayesQuick[[j]],
                          lowerFrequentist = yLowerPoint[[j]],
                          upperFrequentist = yUpperPoint[[j]]
                          )
  }

  class(CI) <- "ci"

  return(CI)
}

plot.ci <- function(CI) {
  #' Plot function for confidence bands
  #'
  #' This function plots the confidence bands
  #' @param CI output from the ci function
  #' @keywords ci print
  #' @export
  J <- length(CI)
  for (j in 1:J) {
    CIpoly <- data.frame(x = c(CI[[j]]$x, rev(CI[[j]]$x)), 
                         y = c(CI[[j]]$lower, rev(CI[[j]]$upper)))

    yRange <- range(c(CI[[j]]$lower, CI[[j]]$upper))

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
