% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcoda.R
\name{gcoda}
\alias{gcoda}
\title{Conditional dependence network inference for compositional data}
\usage{
gcoda(x, counts, pseudo, lambda.min.ratio, nlambda, ebic.gamma)
}
\arguments{
\item{x}{A numeric matrix where rows represent samples and columns represent compositional variables.}

\item{counts}{Logical, indicating whether the compositional data matrix is a count matrix. Default is FALSE.}

\item{pseudo}{Numeric, the pseudo count to add if `counts = TRUE`. Default is 0.5.}

\item{lambda.min.ratio}{Numeric, the ratio of the minimum lambda to the maximum lambda for the lambda sequence. Default is used to create the lambda sequence.}

\item{nlambda}{Integer, the number of lambda values to be used in the lambda sequence. Default is determined by the function.}

\item{ebic.gamma}{Numeric, the EBIC parameter used for model selection. Default is used for EBIC scoring.}
}
\value{
A list containing:
\item{lambda}{The sequence of lambda values used.}
\item{nloglik}{The negative log-likelihood values for each lambda.}
\item{df}{The degrees of freedom for each lambda.}
\item{path}{The list of adjacency matrices for each lambda value.}
\item{icov}{The list of estimated inverse covariance matrices for each lambda value.}
\item{ebic.score}{The EBIC scores for each lambda value.}
\item{opt.index}{The index of the optimal lambda value based on the EBIC score.}
\item{refit}{The refitted adjacency matrix corresponding to the optimal lambda.}
\item{opt.icov}{The optimal estimated inverse covariance matrix.}
\item{opt.lambda}{The optimal lambda value.}
}
\description{
Conditional dependence network inference for compositional data
}
