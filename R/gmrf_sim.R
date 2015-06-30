# -------------------------------------------------------
#
#  Functions to simulate from simple gmrfs
#
# -------------------------------------------------------

#' Simulate from an inhomogenous GMRF
#'
#' Details
#'
#' This function does stuff.
#'
#' @param Q a symmetric positive semi-definate matrix corresponding
#'        to the precision matrix of an inhomogenous GMRF
#' @param tol tolerance (relative to largest eigen value) for numerical lack
#'        of positive-definiteness in `Q'.
#' @param rank the rank of the precision matrix. If this is not supplied it
#'        it is estimated by comparing the eigen values to tol * largest eigenvalue
#' @return a single draw from with the appropriate covariance structure
#' @export
#' @examples
#' ## create a rw2 GMRF precision matrix and simulate
#' Q <- getQrw(100, order = 2)
#' x <- simQ(Q)
#' plot(x)
#' ##
#' ## simulate an AR1 GMRF with toroidal edge correction
#' Q <- getQar(100, phi = 0.8)
#' x <- simQ(exp(5)*Q)
#' plot(x, type = "l")
simQ <- function(Q, tol = 1e-9, rank = NULL) {
  # eigen value decomposition
  eQ <- eigen(Q)
  # simulate
  simeQ(eQ)
}


#' Simulate from an inhomogenous GMRF
#'
#' Required the eigen decomposition of a precision matrix
#'
#' Note if rank is suplied and is less than that true rank, the simulation
#' will be from a reduced rank GMRF.
#'
#' @param eQ an eigen decomposition of a symmetric positive semi-definate matrix corresponding
#'        to the precision matrix of an inhomogenous GMRF
#' @param tol tolerance (relative to largest eigen value) for numerical lack
#'        of positive-definiteness in `Q'.
#' @param rank the rank of the precision matrix. If this is not supplied it
#'        it is estimated by comparing the eigen values to tol * largest eigenvalue.
#' @param k optional scaling of the precision matrix.  This can be used to save
#'        recomputing the eigen decomposition for different smoothing parameters.
#' @return a single draw from with the appropriate conditional covariance structure
#' @export
#' @examples
#' ## create a rw2 GMRF precision matrix and simulate
#' Q <- getQrw(100, order = 2)
#' eQ <- eigen(Q)
#' x <- simeQ(eQ, k = exp(5))
#' plot(x, type = "l")
simeQ <- function(eQ, tol = 1e-9, rank = NULL, k = 1) {
  # of rank not specified, estimate it from the eigen values
  # maybe should check that a suplied rank is sensible
  if (is.null(rank)) rank <- sum(eQ $ values > tol)
  colSums(t(eQ $ vectors[,1:rank]) * rnorm(rank, 0, 1/sqrt(k * eQ $ values[1:rank])))
}
