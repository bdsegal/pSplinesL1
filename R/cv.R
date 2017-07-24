
cv <- function(y,
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
               verbose = 2) {
  #' cross validation for l1-penalized additive mixed models
  #'
  #' This function runs cross-validation one path at a time 
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

  if (is.null(smoothInit)) {
    smoothOpt <- rep(0, J + 1)
  } else {
    smoothOpt <- smoothInit
  }

  names(smoothOpt) <- c("tau", paste("lambda", 1:J, sep = ""))

  # separate warm starts for each fold
  warm <- vector("list", K)

  CVlong <- data.frame(param = factor(rep(c("tau", paste("lambda[", 1:J, "]", sep = "")), 
                                   each = pathLength),
                                     levels = c("tau", 
                                      paste("lambda[", 1:J, "]", sep = ""))),
                       pathVal = do.call(c, smoothPath),
                       cv = NA,
                       sd = NA)

  minInd <- rep(NA, 3)

  if (is.null(paramOrder)) {
    cvOrder <- 1:(J+1)
  } else if(is.numeric(paramOrder) & 
            length(paramOrder) == J + 1 &
            length(unique(paramOrder)) == J + 1) {
    cvOrder <- paramOrder
  } else {
    stop("paramOrder must be either NULL or a vector of length J+1 that 
          contains each of 1 through J+1 exactly once")
  }

  paramNames <- c("tau", paste("lambda", 1:J, sep = ""))

  for (j in cvOrder) {
    for (pathIter in 1:length(smoothPath[[j]])) {
      tauLambda <- smoothOpt
      tauLambda[j] <- smoothPath[[j]][pathIter]
      names(tauLambda) <- c("tau", paste("lambda", 1:J))

      if(verbose >= 1) {
        print("")
        print(paste(paramNames[j], 
                    ", iteration ", pathIter, " of ", 
                    pathLength,  sep = ""))
        print(tauLambda)
      }

      rho <- min(5, max(tauLambda))

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
      CVlong$cv[(j-1)*pathLength + pathIter] <- sum(cvTemp)
      CVlong$sd[(j-1)*pathLength + pathIter] <- sd(cvTemp)

    }
    
    smoothjcv <- CVlong$cv[((j-1)*pathLength + 1):(j*pathLength)]
    smoothjsd <- CVlong$sd[((j-1)*pathLength + 1):(j*pathLength)]
    minInd[j] <- which.min(smoothjcv)
    if (se1) { 
      # TODO: fix
      ind1se <- which.max(smoothjcv <= smoothjcv[minInd[j]] + smoothjsd[minInd[j]])
      smoothOpt[j] <- smoothPath[[j]][ind1se]
    } else {
    # in case of tie, largest parameter selected (smallest df)
      smoothOpt[j] <- smoothPath[[j]][which.min(smoothjcv)]
    }

    cvVlines <- data.frame(pathVal = sapply(1:(J+1), function(l) {
                                           smoothPath[[l]][minInd[l]]}),
                           param = c("tau", paste("lambda[", 1:J, "]", sep = "")))

  }

  gg <- ggplot(aes(x = log(pathVal), y = cv), data = CVlong)+
        geom_point()+
        geom_errorbar(aes(ymax = cv + sd, ymin = cv - sd))+
        facet_wrap(~param, scale = "free", labeller = label_parsed)+
        theme_bw(28)+
        geom_vline(aes(xintercept = log(pathVal)),
                   size = 1, linetype = "dashed", data = cvVlines)+
        labs(x = "log(smoothing parameter)", y = "CV error")

  ret <- list(smoothOpt = smoothOpt, CVlong = CVlong, smoothPath = smoothPath,
              minInd = minInd, cvVlines = cvVlines, gg = gg)
  class(ret) <- "cv"

  return(ret)
}

print.cv <- function(out){
  #' Print function for cv
  #'
  #' This function prints the optimal smoothing parameters
  #' and prints the smoothing parameter paths
  #' @param out output from the cv function
  #' @keywords cv print
  #' @export
  
  print(out$smoothOpt)
  print(out$gg)
}
