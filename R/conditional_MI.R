#' Conditional Mutual Information Network Construction
#'
#' Constructs a network using conditional mutual information (CMI) to identify relationships between taxa.
#' The function applies two levels of CMI computation, using quantile thresholds to filter edges in the resulting adjacency matrices.
#'
#' @param data A numeric matrix where rows represent taxa and columns represent samples. If \code{quantitative} is TRUE, data will be log-transformed.
#' @param q1 A numeric value representing the quantile threshold for filtering edges in the order 0 adjacency matrix.
#' @param q2 A numeric value representing the quantile threshold for filtering edges in the order 1 adjacency matrix.
#' @param quantitative A logical value indicating if the data is quantitative. If TRUE, data is log-transformed.
#'
#' @return A list containing:
#'   \item{G_order0}{The adjacency matrix after order 0 calculation.}
#'   \item{G_order1}{The adjacency matrix after order 1 calculation.}
#'   \item{Gval_order0}{The CMI values for each taxa pair in the order 0 adjacency matrix.}
#'   \item{Gval_order1}{The CMI values for each taxa pair in the order 1 adjacency matrix.}
#'   \item{quantile_order0}{The quantile threshold used for filtering the order 0 adjacency matrix.}
#'   \item{quantile_order1}{The quantile threshold used for filtering the order 1 adjacency matrix.}
#'   \item{sum_order0}{The sum of edges in the order 0 adjacency matrix.}
#'   \item{sum_order1}{The sum of edges in the order 1 adjacency matrix.}
#'
#' @details The function computes the conditional mutual information for taxa pairs, applying two levels of calculation:
#'   \enumerate{
#'     \item \strong{Order 0}: Computes CMI between each pair of taxa and filters edges based on \code{q1}.
#'     \item \strong{Order 1}: Refines the network by conditioning on shared neighbors, filtering edges based on \code{q2}.
#'   }
#' @importFrom stats cov quantile var
#' @export

conditional_MI <- function(data, q1, q2, quantitative) {
  # Transform data based on quantitative input
  if (quantitative) {
    data <- t(log(data + 1))
  } else {
    data <- t(data)
  }

  n_taxa <- nrow(data)
  G <- matrix(1, n_taxa, n_taxa)
  diag(G) <- 0
  Gval <- G

  # Conditional Mutual Information function
  cmi <- function(v1, v2, vcs = NULL) {
    if (is.null(vcs)) {
      c1 <- var(v1)
      c2 <- var(v2)
      c3 <- det(cov(cbind(v1, v2)))
      cmiv <- 0.5 * log((c1 * c2) / c3)
    } else {
      c1 <- det(cov(cbind(v1, vcs)))
      c2 <- det(cov(cbind(v2, vcs)))
      c3 <- det(cov(cbind(vcs)))
      c4 <- det(cov(cbind(v1, v2, vcs)))
      cmiv <- 0.5 * log((c1 * c2) / (c3 * c4))
    }

    if (!is.finite(cmiv)) {
      cmiv <- 0
    }
    return(cmiv)
  }

  # Order 0 calculation
  for (i in 1:(n_taxa - 1)) {
    for (j in (i + 1):n_taxa) {
      if (!is.na(G[i, j]) && G[i, j] != 0) {

        Gval[i, j] <- cmi(data[i, ], data[j, ])
        Gval[j, i] <- Gval[i, j]
      }
    }
  }

  quantile_order0 <- quantile(Gval, probs = q1)
  G[Gval < quantile_order0] <- 0
  G_order0 <- G
  Gval_order0 <- Gval

  # Order 1 calculation
  for (i in 1:(n_taxa - 1)) {
    for (j in (i + 1):n_taxa) {
      if (!is.na(G[i, j]) && G[i, j] != 0) {

        adj <- which(G[i, ] != 0 & G[j, ] != 0)
        if (length(adj) >= 1) {
          cmiv <- 0
          v1 <- data[i, ]
          v2 <- data[j, ]
          for (k in 1:length(adj)) {
            vcs <- data[adj[k], ]
            cmiv <- max(cmiv, cmi(v1, v2, vcs))
          }
          Gval[i, j] <- cmiv
          Gval[j, i] <- cmiv
        }
      }
    }
  }

  quantile_order1 <- quantile(Gval, probs = q2)
  G[Gval < quantile_order1] <- 0
  G_order1 <- G
  Gval_order1 <- Gval

  # Output
  out <- list(
    G_order0 = G_order0,
    G_order1 = G_order1,
    Gval_order0 = Gval_order0,
    Gval_order1 = Gval_order1,
    quantile_order0 = quantile_order0,
    quantile_order1 = quantile_order1,
    sum_order0 = sum(G_order0),
    sum_order1 = sum(G_order1)
  )
  return(out)
}
