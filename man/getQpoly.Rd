% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gmrf_Qfuncs.R
\name{getQpoly}
\alias{getQpoly}
\title{Differences defining a regional model from spatial polygons...}
\usage{
getQpoly(poly, weights = NULL)
}
\arguments{
\item{poly}{a spatial polygon}

\item{weights}{weights to be applied to the node differences in effect
allowing different connections between regions to vary differently.
An example of a weight could be the length of the shared border between regions}
}
\value{
a Matrix with a column for each region and a row for each connection / edge
}
\description{
Details
}
\details{
Description - This function does stuff.
}

