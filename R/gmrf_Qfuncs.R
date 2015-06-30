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
    if (length(weights) != n - 1 + cyclic) {
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
#' @export
getQrw <- function(n, order = 2, weights = NULL, cyclic = FALSE) {
  D <- Drw(n, order = order, cyclic = cyclic)
  if (!is.null(weights)) {
    if (length(weights) != n - order + order*cyclic) {
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
    if (length(weights) != n) {
      stop("weights must be the same length as the number of unique differences in the GMRF")
    }
    out <- t(D) %*% diag(weights) %*% D
  } else {
    out <- t(D) %*% D
  }
  out
}