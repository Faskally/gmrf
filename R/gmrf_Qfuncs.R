# -------------------------------------------------------
#
#  Functions to create simple gmrfs
#
# -------------------------------------------------------


# -------------------------------
# first the difference functions
# -------------------------------

#' Differences defining a cyclic auto regressive process
#' or equivalently, an AR with toroidal edge effects.
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param phi the autoregressive parameters
#' @return a Matrix
#' @export
Dar <- function(n, phi) {
  order <- length(phi)
  vec <- c(1, rep(0, n-order-1), -rev(phi))
  mat <- sapply(1:n, function(i) vec[(-i + 1 + 0:(n-1)) %% n + 1])
  mat
}


#' Differences defining a 1st order random walk
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF.
#' @param cyclic logical: should the differences be cyclic, i.e.
#'        corrected for toroidal edge effects .
#' @return a Matrix
#' @export
Drw1 <- function(n, cyclic = FALSE) {
    out <- diag(n) 
    diag(out[,-1]) <- -1
    out[n,1] <- -1
    if (cyclic) out else out[-n,]
}

#' Differences defining a cyclic nth order random walk
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF.
#' @param order the order of the random walk.
#' @param cyclic logical: should the differences be cyclic, i.e.
#'        corrected for toroidal edge effects. 
#' @return a Matrix
#' @export
Drw <- function(n, order = 2, cyclic = FALSE) {
  stopifnot(order > 0)
  out <- diag(n)
  for (i in 1:order) {
    out <- out %*% Drw1(n, cyclic = TRUE)
  }
  if (cyclic) out else out[1:(n-order),]
}


#' Differences defining a cyclic nth order random walk
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param wavelength the wavelength of the harmonic. default is for this
#'        this to be n, i.e. one cycle.
#' @param cyclic logical: should the smoother be cyclic, i.e.
#'        corrected for toroidal edge effects.      
#' @return a Matrix
#' @export
Dharm <- function(n, wavelength = n, cyclic = FALSE) {
  D1 <- Drw1(n, cyclic = TRUE)
  D3 <- D1 %*% D1 %*% D1
  out <- D3 + (2*pi/wavelength)^2 * D1
  if (cyclic) out else out[1:(n-3),]
}


#' Differences defining a seasonal model
#'
#' This model treats the sum over m time points
#' to be stationary
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param m the length of the seasonal period
#' @return a Matrix
#' @export
Dseasonal <- function(n, m) {
  vec <- c(rep(1, m), rep(0, n-m))
  t(sapply(1:(n-m+1), function(i) vec[(-i + 1 + 0:(n-1)) %% n + 1]))
}


#' Differences defining a regional model
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb the neighbourhood structure of a regional layout
#' @return a Matrix with a column for each region and a row for each connection / edge
#' @export
Dnb <- function(nb) {
  nnodes <- length(nb)
  # keep edges in one direction only
  nb <- lapply(1:nnodes, function(i) nb[[i]] <- nb[[i]][nb[[i]] > i])
  edge_to <- unlist(nb)
  nedges <- length(edge_to)
  nedge_from <- sapply(nb, length)
  out <- diag(nnodes)[rep(1:nnodes, nedge_from),]
  ids <- cbind(1:nedges, unlist(nb))
  out[ids] <- -1
  out
}



# -------------------------------
# Now the precision functions
# -------------------------------


#' Compute RW1 precision matrix
#'
#' If weights are supplied
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param weights weights to be applied to the node differences (see details)
#' @param cyclic logical: should the smoother be cyclic, i.e.
#'        corrected for toroidal edge effects. 
#' @return what does it return
#' @export
getQrw1 <- function(n, weights = NULL, cyclic = FALSE) {
  D <- Drw1(n, cyclic = cyclic)
  if (!is.null(weights)) {
    if (length(weights) != nrow(D)) {
      stop("weights must be the same length as the number of unique differences in the GMRF")
    }
    out <- t(D) %*% diag(weights) %*% D
  } else {
    out <- t(D) %*% D
  }
  out
}




#' Compute RWn precision matrix
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param order the order of the random walk 
#' @param weights weights to be applied to the node differences (see details)
#' @param cyclic logical: should the smoother be cyclic, i.e.
#'        corrected for toroidal edge effects. 
#' @return what does it return
#' @examples
#' require(gmrf)
#' n <- 100
#' idx <- 1:n
#' idy <- c(5:40, 60:n)
#' set.seed(64)
#' Q <- getQrw(n, order = 2)
#' x <- simQ(exp(1) * Q)
#' # simulate the first 3/4, say
#' y <- x + rnorm(n) * 3.5
#' y <- y[idy]
#' ## set up variables for smoothing
#' rownames(Q) <- colnames(Q) <- idx
#' ## fit an RW2 smoother with restricted df
#' g1 <- gam(y ~ s(idy, bs = "gmrf", xt = list(penalty = Q), k = length(y)-1), method="REML")
#' summary(g1)
#' plot(idy, y, xlim = range(idx), ylim = range(x,y))
#' lines(id, x, col = "blue", lwd = 2)
#' pred <- predict(g1, newdata = list(idy = idx), se = TRUE)
#' lines(idx, pred$fit, col = "red", lwd = 2)
#' lines(idx, pred$fit + 2*pred$se.fit, col = "red", lty = 2)
#' lines(idx, pred$fit - 2*pred$se.fit, col = "red", lty = 2)
#' @export
getQrw <- function(n, order = 2, weights = NULL, cyclic = FALSE) {
  D <- Drw(n, order = order, cyclic = cyclic)
  if (!is.null(weights)) {
    if (length(weights) != nrow(D)) {
      stop("weights must be the same length as the number of unique differences in the GMRF")
    }
    out <- t(D) %*% diag(weights) %*% D
  } else {
    out <- t(D) %*% D
  }
  out
}


#' Differences defining a cyclic auto regressive process
#' or equivalently, an AR with toroidal edge effects.
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param n the size of the GMRF
#' @param phi the autoregressive parameters
#' @param weights weights to be applied to the node differences in effect
#'        allowing different variances at each time step.
#' @return a Matrix
#' @export
getQar <- function(n, phi, weights = NULL) {
  D <- Dar(n, phi)
  if (!is.null(weights)) {
    if (length(weights) != nrow(D)) {
      stop("weights must be the same length as the number of unique differences in the GMRF")
    }
    out <- t(D) %*% diag(weights) %*% D
  } else {
    out <- t(D) %*% D
  }
  out
}


#' Differences defining a regional model
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb the neighbourhood structure of a regional layout
#' @param weights weights to be applied to the node differences in effect
#'        allowing different connections between regions to vary differently.
#'        An example of a weight could be the length of the shared border between regions
#' @return a Matrix with a column for each region and a row for each connection / edge
#' @export
getQnb <- function(nb, weights = NULL) {
  D <- Dnb(nb)
  if (!is.null(weights)) {
    if (length(weights) != nrow(D)) {
      stop("weights must be the same length as the number of unique differences in the GMRF")
    }
    out <- t(D) %*% diag(weights) %*% D
  } else {
    out <- t(D) %*% D
  }
  out
}
