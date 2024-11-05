################################################################################
# File: SparCC.R
# Aim : SparCC
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 11/12/2014
#-------------------------------------------------------------------------------
# SparCC for fractions known
#   function: SparCC.frac
#   input:
#          x ------ nxp fraction data matrix, row is sample, col is variable
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#' SparCC for Fraction Data
#'
#' This function performs SparCC analysis for fraction data to estimate covariance and correlation matrices.
#'
#' @param x A numeric matrix where rows represent samples and columns represent variables (fraction data).
#' @param kmax Integer, the maximum iteration steps for SparCC. Default is 10.
#' @param alpha Numeric, the threshold for strong correlation. Default is 0.1.
#' @param Vmin Numeric, the minimal variance if negative variance appears. Default is 1e-4.
#'
#' @return A list containing:
#' \item{cov.w}{The estimated covariance matrix.}
#' \item{cor.w}{The estimated correlation matrix.}
#' @importFrom stats median

#' @export
SparCC.frac <- function(x, kmax, alpha, Vmin) {
  # Log transformation
  x <- log(x);
  p <- ncol(x);
  # T0 = var(log(xi/xj)) variation matrix
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  #cor.w <- (outer(var.w, var.w, "+") - T0 ) /
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  Is <- sqrt(1/var.w);
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor.w[cor.w <= - 1] <- - 1;
  cor.w[cor.w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1;
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor.w <- cor.w;
    # Remove diagonal
    diag(curr_cor.w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor.w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) /
      #  sqrt(outer(var.w, var.w, "*")) / 2;
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) *
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1;
      cor.w[cor.w >= 1] <- 1;
    }
    else {
      break;
    }
    #
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov.w <- cor.w * Is * rep(Is, each = p);
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}
#-------------------------------------------------------------------------------
