% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SparCC.frac.R
\name{SparCC.frac}
\alias{SparCC.frac}
\title{SparCC for Fraction Data}
\usage{
SparCC.frac(x, kmax, alpha, Vmin)
}
\arguments{
\item{x}{A numeric matrix where rows represent samples and columns represent variables (fraction data).}

\item{kmax}{Integer, the maximum iteration steps for SparCC. Default is 10.}

\item{alpha}{Numeric, the threshold for strong correlation. Default is 0.1.}

\item{Vmin}{Numeric, the minimal variance if negative variance appears. Default is 1e-4.}
}
\value{
A list containing:
\item{cov.w}{The estimated covariance matrix.}
\item{cor.w}{The estimated correlation matrix.}
}
\description{
This function performs SparCC analysis for fraction data to estimate covariance and correlation matrices.
}
