################################################################################
# File: SparCC.R
# Aim : SparCC
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 11/12/2014
#-------------------------------------------------------------------------------
# SparCC for counts known
#   function: SparCC.count
#   input:
#          x ------ nxp count data matrix, row is sample, col is variable
#       imax ------ resampling times from posterior distribution. default 20
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#' SparCC for Counts Data
#'
#' This function performs SparCC analysis for count data to estimate covariance and correlation matrices.
#'
#' @param x A numeric matrix where rows represent samples and columns represent variables (count data).
#' @param imax Integer, the number of resampling times from the posterior distribution. Default is 20.
#' @param kmax Integer, the maximum iteration steps for SparCC. Default is 10.
#' @param alpha Numeric, the threshold for strong correlation. Default is 0.1.
#' @param Vmin Numeric, the minimal variance if negative variance appears. Default is 1e-4.
#'
#' @return A list containing:
#' \item{cov.w}{The estimated covariance matrix.}
#' \item{cor.w}{The estimated correlation matrix.}
#' @importFrom stats median
#' @importFrom gtools rdirichlet

#' @export

SparCC.count <- function(x, imax, kmax, alpha, Vmin) {
  # dimension for w (latent variables)
  p <- ncol(x);
  n <- nrow(x);
  # posterior distribution (alpha)
  x <- x + 1;
  # store generate data
  y <- matrix(0, n, p);
  # store covariance/correlation matrix
  cov.w <- cor.w <- matrix(0, p, p);
  indLow <- lower.tri(cov.w, diag = T);
  # store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x)
      gtools::rdirichlet(n = 1, alpha = x)));
    # estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # store variance/correlation only low triangle
    covs[, i] <- cov_cor$cov.w[indLow];
    cors[, i] <- cov_cor$cor.w[indLow];
  }
  # calculate median for several posterior samples
  cov.w[indLow] <- apply(covs, 1, median);
  cor.w[indLow] <- apply(cors, 1, median);
  #
  cov.w <- cov.w + t(cov.w);
  diag(cov.w) <- diag(cov.w) / 2;
  cor.w <- cor.w + t(cor.w);
  diag(cor.w) <- 1;
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}
