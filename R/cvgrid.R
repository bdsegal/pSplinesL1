
cvgrid <- function(y,
               X,
               rand,
               id,
               data,
               K = 5,
               se1 = FALSE,
               smoothInit = NULL,
               pathLength = 20,
               useTestRanEff = TRUE,
               tauMax = NULL,
               revPath = FALSE,
               paramOrder = NULL,
               epsilonAbs = 1e-4,
               epsilonRel = 1e-4,
               iterMax = 1e3,
               centerZ = FALSE,
               verbose = 2,
               fileName = NULL) {
  #' cross validation grid search for l1-penalized additive mixed models
  #'
  #' This function runs cross-validation grid search
  #' to estimate the smoothing parameters
  #' @param y character string denoting output variable
  #' @param X list of fixed effect design matrices setup by l1smooth
  #' @param rand list of random effect terms setup by randomEffects
  #' @param id character string denoting subject id variable
  #' @param data dataset
  #' @param se1 TRUE/FALSE if true, sets smoothing parameter to value with
  #'        maximum CV error within 1 sd of the value that minimizes CV error
  #' @param smoothInit vector of initial smoothing parameter values
  #' @param pathLength length of each smoothing parameter path
  #'        (evenly spaced on the log scale)
  #' @param useTestRanEff TRUE/FALSE if true, estimates random effects with test sample
  #' @param tauMax scalar max value to use in tau path. If NULL, set to maximum
  #'        value used for fixed effect smoothing parameters.
  #' @param revPath TRUE/FALSE if true, estimates smoothing paths from smallest to largest
  #' @param paramOrder order of smoothing parameters to be fit
  #' @param epsilonAbs absolute error in ADMM (tolerance for convergence)
  #' @param epsilonRel relative error in ADMM (tolerance for convergence)
  #' @param iterMax maximum number of ADMM iterations
  #' @param centerZ TRUE/FALSE If True, imposes centering constraints on random effects
  #' @param verbose integer (0, 1, 2) amount of information printed to terminal. 0 = no information, 1 = path iteration, 2 = path iteration + fold results
  #' @param saveFile file name for saving CV results as csv file, e.g. "CV_results". Do not include extension. If NULL, no file saved.
  #' @keywords ADMM, cross validation
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
  #' # run cross validation and view results
  #' cvOut <- cvgrid(y = "y", id = "id", X = X, rand = rand,
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

  ycol <- which(colnames(data) == y)
  J <- length(X) # number of smooths

  lambdaMax <- rep(NA, J)
  for (j in 1:J) {
    lambdaMax[j] <- max(abs(solve(tcrossprod(X[[j]]$D, X[[j]]$D), 
                             tcrossprod(X[[j]]$D, X[[j]]$F) %*% data[, ycol])))
  }

  # By default use largest lambdaMax for tau (random effects)
  # This might not be appropriate for all data
  if(is.null(tauMax)) {
    tauLambdaMax <- c(max(lambdaMax), lambdaMax)
  } else {
      tauLambdaMax <- c(tauMax, lambdaMax)
  }

  # TODO: allow for different smoothing paths values
  smoothPath <- lapply(tauLambdaMax, function(x){
                       exp(seq(log(x), 
                               log(x*1e-5), 
                       length.out = pathLength))})
  if (revPath) {
    smoothPath <- lapply(smoothPath, rev)
  }
  
  idCol <- which(colnames(data) == id)
  uniqId <- unique(data[, idCol])
  
  # if indicator variable, setup folds with subjects in level
  byVars <- unlist(sapply(X, function(x) {x$by}))
  if (is.null(byVars)) {
    # Create folds without regard to having subjects from each cell in each fold
    fold <- list()
    nid <- floor(length(uniqId) / K)
    for (k in 1:K) {
      if (k < K) {
        fold[[k]] <- uniqId[((k-1)*nid + 1):(k*nid)]
      } else {
        fold[[k]] <- uniqId[((k-1)*nid + 1):length(uniqId)]
      }
    }
  } else {

    # split cells evenly across folds
    byCols <- which(colnames(data) %in% byVars)
    byIsInd <- sapply(byCols, function(i) {
                         length(unique(data[, i])) == 2})
    
    # get variable patterns for each id
    varPatterns <- matrix(NA, nrow = length(uniqId), ncol = sum(byIsInd))
    rownames(varPatterns) <- uniqId
    colnames(varPatterns) <- byVars[byIsInd]
    for(i in 1:length(uniqId)) {
      varPatterns[i, ] <- data[which(data[, idCol] == uniqId[i])[1], 
                       c(byCols[byIsInd])]
    }

    # get unique variable patterns
    uniqVarPatterns <- as.data.frame(unique(varPatterns))
    rownames(uniqVarPatterns) <- NULL

    # get ids with each unique variable pattern
    ids <- list()
    for(i in 1:nrow(uniqVarPatterns)) {
      ind <- rep(TRUE, nrow(varPatterns))
      for (d in 1:ncol(uniqVarPatterns)) {
        ind <- ind * (varPatterns[, d] == uniqVarPatterns[i, d])
      }
      ids[[i]] <- rownames(varPatterns)[which(ind == 1)]
    }

    # to do: give warning if not enough subjects in each cell
    minCellInd <- which.min(sapply(ids, length))
    minCell <- min(sapply(ids, length))
    if(floor(minCell / K) < 2) {
      stop(paste("Not enough subjects for each combination of factor levels 
         to have at least 2 in each fold (only ", minCell, " subjects in smallest
         cell). Try reducing the number of folds to K <= ", floor(minCell/2),".", sep = "")
          )
    }
    # tab <- xtabs(as.formula(paste("~", paste(colnames(varPatterns), collapse = " + "), sep = " ")), data = mat)

    fold <- list()
    nid <- floor(sapply(ids, length) / K)
    for (k in 1:K) {
      if (k < K) {
        fold[[k]] <- do.call(c, 
                  sapply(1:length(ids), function(i) {
                    ids[[i]][((k-1)*nid[i] + 1):(k*nid[i])]
                    }))
      } else {
        fold[[k]] <- do.call(c, 
                  sapply(1:length(ids), function(i) {
                    ids[[i]][((k-1)*nid[i] + 1):length(ids[[i]])]
                    }))
      }
    }
  }

  # separate warm starts for each fold
  warm <- vector("list", K)

  CVlong <- expand.grid(smoothPath)
  colnames(CVlong) <- c("tau", paste("lambda[", 1:J, "]", sep = ""))
  CVlong$cv <- NA
  CVlong$sd <- NA

  fileName <- paste(fileName, ".csv", sep = "")

  for (gridIter in 1:nrow(CVlong)) {

    if(verbose >= 1) {
      print("")
      print(paste("iteration ", gridIter, " of ", 
                  nrow(CVlong),  sep = ""))
      print(CVlong[gridIter, 1:(J+1)])
    }

    rho <- min(5, max(CVlong[gridIter, 1:(J+1)]))

    # k-fold cross-validation
    cvTemp <- rep(0, K)
    for (k in 1:K) {
      testRec <- (data$id %in% fold[[k]])
      train <- data[!testRec, ]
      test <- data[testRec, ]

      nTrain <- table(train$id)
      nTest <- table(test$id)

      Xtrain <- Xtest <- X
      for (l in 1:J) {
        Xtrain[[l]]$F <- X[[l]]$F[!testRec, ]
        Xtest[[l]]$F <- X[[l]]$F[testRec, ]
      }

      Ztrain <- rand$Z[which(!rownames(rand$Z) %in% fold[[k]]), 
                  which(!colnames(rand$Z) %in% fold[[k]])]

      Strain <- rand$S[which(!rownames(rand$S) %in% fold[[k]]), 
                  which(!colnames(rand$S) %in% fold[[k]])]           

      tauLambda <- as.numeric(CVlong[gridIter, 1:(J+1)])
      m1 <- admm(y = y, X = Xtrain, Z = Ztrain, S = Strain,
           lambda = tauLambda[2:(J + 1)],
           tau = tauLambda[1],
           rho = rho,
           epsilonAbs = epsilonAbs,
           epsilonRel = epsilonRel,
           iterMax = iterMax,
           warm = warm[[k]],
           data = train,
           centerZ = centerZ,
           forCV = TRUE)

      yTestHat <- m1$warm$beta0
      for (l in 1:J) {
        yTestHat <- yTestHat + Xtest[[l]]$F %*% m1$warm$beta[[l]]
      }

      if (useTestRanEff) {
        Ztest <- rand$Z[which(rownames(rand$Z) %in% fold[[k]]), 
                  which(colnames(rand$Z) %in% fold[[k]])]

        Stest <- rand$S[which(rownames(rand$S) %in% fold[[k]]), 
                  which(colnames(rand$S) %in% fold[[k]])] 
        bTest <- solve(crossprod(Ztest) + tauLambda[1] * Stest,
                       crossprod(Ztest, test[, ycol] - yTestHat))
        yTestHat <- yTestHat + Ztest %*% bTest
      }

      cvTemp[k] <- sum((test[, ycol] - yTestHat)^2)

      warm[[k]] <- m1$warm

      if(verbose >= 2) {
        print(paste("fold: ", k,
                ", ADMM iterations: ", m1$conv$iter, 
                ", cumulative CV error: ", signif(sum(cvTemp), 4), 
                ", training error: ", signif(sum(m1$fit$residuals^2), 4),
                sep = ""))
      }
    }

    CVlong$cv[gridIter] <- sum(cvTemp)
    CVlong$sd[gridIter] <- sd(cvTemp)

    # write CV results to file on iteration at a time
    if(!is.null(fileName)) {
      if(gridIter == 1) {
        write.table(CVlong[1, ], file = fileName, row.names = FALSE,
                    sep = ",")
      } else {
        write.table(CVlong[gridIter, ], file = fileName, row.names = FALSE,
                    sep = ",", col.names = FALSE, append = TRUE)
      }
    }

  }
    
  minInd <- which.min(CVlong$cv)
  if (se1) { 
    ind1se <- which.max(CVlong$cv <= CVlong$cv[minInd] + CVlong$sd[minInd])
    smoothOpt <- as.numeric(CVlong[ind1se, 1:(J+1)])
  } else {
    # in case of tie, largest parameter selected (smallest df)
    smoothOpt <- as.numeric(CVlong[minInd, 1:(J+1)])
  }


  ret <- list(smoothOpt = smoothOpt, CVlong = CVlong, smoothPath = smoothPath,
              minInd = minInd)
  class(ret) <- "cvgrid"

  return(ret)
}

print.cvgrid <- function(out){
  #' Print function for cvgrid
  #'
  #' This function prints the optimal smoothing parameters
  #' and prints the smoothing parameter paths
  #' @param out output from the cvgrid function
  #' @keywords cvgrid print
  #' @export
  
  print(out$smoothOpt)
}
