
# function for setting up splines --------------------------------------------_
ps <- function(x, norder, k, data, by = NULL, ref = NULL, width = NULL, delta = 1e-7,
                     center = TRUE) {
  #' Function for setting up p-spline smooths
  #'
  #' @param x vector evaluation points (typically time)
  #' @param norder scalar order of B-spline bases (order k is degree k-1)
  #' @param k scalar k+1 is order of finite difference matrix
  #' @param data dataset in which to evaluate arguments
  #' @param by optional character string designating varying coefficient
  #' @param ref optional reference level of by variable 
  #' @param width scalar width of B-spline bases
  #' @param delta scalar place upper knot in B-spline basis at max(x) + delta
  #' @param center TRUE/FALSE if True, impose centering constrain on non-variable coefficient models
  #' @keywords additive components
  #' @export
  #' @examples
  #' library(psplinesl1)
  #' data(simData)
  #'
  #' # setup p-spline matrices
  #' X <- list(ps(x = "x", data = simData, 
  #'              norder = 2, k = 1, width = 0.05,
  #'              center = TRUE))

  xcol <- which(colnames(data) == x)

  xRange <- c(floor(min(data[, xcol])), ceiling(max(data[, xcol])) + delta)
  if(is.null(width)) width <- diff(xRange) / 10
  basis <- create.bspline.basis(rangeval = xRange,
                                breaks = seq(xRange[1], xRange[2], width),
                                norder = norder)
  F <- eval.basis(data[, xcol], basis)
  p <- basis$nbasis
  D <- diff(diag(p), diff = k + 1)

  if (is.null(by)) {
    if(center) {
      # centering constraints
      Ct <- as.matrix(colSums(F), nrow = ncol(F))
      Q <- qr.Q(qr(Ct), complete = TRUE)[, -1]
      Fq <- F %*% Q
      Dq <- D %*% Q

      return(list(basis = basis, F = Fq, Q = Q, p = p,
                  D = Dq, by = by, k = k, width = width, center = center))
    } else {

      return(list(basis = basis, F = F, Q = NULL, p = p,
                  D = D, by = by, ref = ref, k = k, width = width, center = center))
    }
  }

  else if (!is.null(by)) {

    byCol <- which(colnames(data) == by)
    if(length(byCol) == 0) {
      stop(paste("The 'by' variable '", by, "' not found in the data.", sep = ""))
    } else if(length(byCol) > 1) {
      stop(paste("The 'by' variable '", by, "' matches ", length(byCol),
                 " variables in the data.", sep = ""))
    } else {

      if(is.factor(data[, byCol])) {

        if(length(unique(data[, byCol])) != 2 ) {
          stop("The 'by' argument only supports factors with two levels, as well
         as numeric variables. To include a factor with >2 levels, please
         make separate numeric (0/1) indicator variables for all levels
         except the reference level, and make separate l1smooth objects 
         for each indicator variable.")
        }

        if(is.null(ref)) {
          ref <- levels(data[, byCol])[1]
        }
        refVar <- paste(by, ref, sep = "")
        byMat <- model.matrix(~ -1 + type, data = data)
        byVec <- byMat[, -which(colnames(byMat) == refVar)]
      } else {
        byVec <- data[, byCol]
      }

      if(is.factor(length(unique(data[, byCol])) < 2)) {
        stop("'by' variable must have >1 unique value")
      }
      F <- sweep(F, 1, byVec, `*`)

      return(list(basis = basis, F = F, Q = NULL, p = p,
                  D = D, by = by, ref = ref, k = k, width = width, center))
    }
  }
}

re <- function(x, id, data, randomCurves = FALSE, width = 5, 
                          norder = 4, derivOrder = 2) {
  #' Function for setting up random effects (intercepts or B-spline curves)
  #'
  #' @param x vector evaluation points (typically time)
  #' @param id character string designating grouping variable
  #' @param data dataset in which to evaluate arguments
  #' @param randomCurves TRUE/FALSE if True, uses B-spline bases, otherwise random intercepts
  #' @param norder scalar order of B-spline bases (order k is degree k-1) (if randomCurves = TRUE)
  #' @param width scalar width of B-spline bases (if randomCurves = TRUE)
  #' @param derivOrder scalar order of derivative used in the random effect penalty matrix (if randomCurves = TRUE)
  #' @keywords random effects
  #' @export
  #' @examples
  #' library(psplinesl1)
  #' data(simData)
  #'
  #' # setup random effect matrices
  #' rand <- re(x = "x", id = "id", data = simData,
  #'            randomCurves = FALSE)
  
  n <- table(data$id)
  idcol <- which(colnames(data) == id)
  Zlist <- list()

  if(randomCurves) {
      
    xcol <- which(colnames(data) == x)

    xRange <- c(floor(min(data[, xcol])), ceiling(max(data[, xcol])))
    basisZ <- create.bspline.basis(rangeval = xRange,
                                    breaks = seq(xRange[1], xRange[2], width),
                                    norder = norder)
    Zlong <- eval.basis(data[, xcol], basisZ)
    cumn <- c(0,cumsum(n))
    for (i in 1:length(n)) {
      Zlist[[i]] <- Zlong[(cumn[i]+1):cumn[i+1], ]
    }
    Z <- bdiag(Zlist)

    Stemp <- bsplinepen(basisZ, Lfdobj = derivOrder)
    Slist <- list()
    for (i in 1:length(n)) {
      Slist[[i]] <- Stemp
    }
    S <- bdiag(Slist)

    colnames(Z) <- rep(unique(data[, idcol]), times = sapply(Zlist, ncol))
    rownames(Z) <- data[, idcol]

    colnames(S) <- rownames(S) <- 
      rep(unique(data[, idcol]), times = sapply(Slist, ncol))
  } else {

    for (i in 1:length(n)) {
      Zlist[[i]] <- rep(1, n[i])
    }
    Z <- as.matrix(bdiag(Zlist))
    S <- diag(ncol(Z))

    rownames(Z) <- data[, idcol]
    colnames(Z) <- unique(data[, idcol])
    colnames(S) <- rownames(S) <- unique(data[, idcol])
  }

  return(list(Z = as.matrix(Z), S = as.matrix(S)))

}