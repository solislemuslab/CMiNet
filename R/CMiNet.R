utils::globalVariables(c("enabled", "params", "Pair", "Hamming_Distance", "Common_Edges", "Edges_Net1", "Edges_Net2", "Weight", "name"))
#' This function constructs a consensus network from microbiome data using multiple methods.
#'
#' @param data A numeric matrix or data frame of microbiome data where rows represent samples and columns represent features.
#' @param quantitative A logical value indicating if the data is quantitative. Defaults to TRUE.
#' @param TT A numeric value representing the quantile threshold used to binarize the adjacency matrices.
#' @param pearson A list with parameters for Pearson correlation (`enabled`, `params`).
#' @param spearman A list with parameters for Spearman correlation (`enabled`, `params`).
#' @param bicor A list with parameters for biweight midcorrelation (`enabled`, `params`).
#' @param sparcc A list with parameters for SparCC correlation (`enabled`, `params`).
#' @param spiecEasi_mb A list with parameters for SpiecEasi Meinshausen-Bühlmann method (enabled, params).
#' @param spiecEasi_glasso A list with parameters for SpiecEasi Graphical Lasso method (enabled, params).
#' @param spring A list with parameters for SPRING method (`enabled`, `params`).
#' @param gcoda A list with parameters for GCODA method (`enabled`, `params`).
#' @param c_MI A list with parameters for Conditional Mutual Information method (`enabled`, `params`).
#' @param cclasso A list with parameters for CCLasso method (`enabled`, `params`).
#' @importFrom SpiecEasi spiec.easi
#' @return A list containing the consensus network, edge list, and any errors encountered during network construction.
#' @importFrom WGCNA bicor
#' @importFrom SPRING SPRING
#' @importFrom stats cor quantile var pt rnorm
#' @import readr
#' @importFrom utils write.csv
#' @export
CMiNet <- function(data, quantitative = TRUE, TT,
                                         pearson = list(enabled , params = list()),
                                         spearman = list(enabled, params = list()),
                                         bicor = list(enabled, params = list()),
                                         sparcc = list(enabled, params),
                                         spiecEasi_mb = list(enabled, params),
                                         spiecEasi_glasso = list(enabled, params),
                                         spring = list(enabled, params ),
                                         gcoda = list(enabled, params ),
                                         c_MI = list(enabled, params),
                                         cclasso = list(enabled, params)) {

  # Convert data to matrix to avoid type issues
  data <- as.matrix(data)

  # Create folders named 'Network' and 'Binary_Network' to save the matrices
  if (!dir.exists("Network")) {
    dir.create("Network")
  }
  if (!dir.exists("Binary_Network")) {
    dir.create("Binary_Network")
  }

  # Create a file to save errors and warnings
  error_log <- file("Network/errors_warnings.txt", open = "wt")

  networks <- list()
  errors <- list()
  num_features <- ncol(data)

  construct_network <- function(method, adj_matrix, use_quantile = TRUE) {
    # Save the adjacency matrix before binarization
    write.csv(adj_matrix, file.path("Network", paste0(method, ".csv")), row.names = TRUE)

    if (use_quantile) {
      # Threshold the matrix based on quantile TT
      th <- quantile(abs(adj_matrix), probs = TT, na.rm = TRUE)
      binary_matrix <- ifelse(abs(adj_matrix) >= th, 1, 0)
    } else {
      # Set non-zero values to 1 and others to 0
      binary_matrix <- ifelse(adj_matrix != 0, 1, 0)
    }
    # Save the binary matrix
    write.csv(binary_matrix, file.path("Binary_Network", paste0(method, "_binary.csv")), row.names = TRUE)
    return(list(adj_matrix = adj_matrix, binary_matrix = binary_matrix))
  }

  # Suppress warnings and messages
  suppressWarnings(suppressMessages(
    {
      if (pearson$enabled) {
        start_time <- Sys.time()
        tryCatch({
          adj_matrix <- cor(data, method = "pearson")
          if (is.null(adj_matrix)) stop("Pearson returned NULL")
          networks$pearson <- construct_network("pearson", adj_matrix, use_quantile = TRUE)
        }, error = function(e) {
          errors$pearson <- e$message
          writeLines(paste("Error in Pearson:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing Pearson correlation network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (spearman$enabled) {
        start_time <- Sys.time()
        tryCatch({
          adj_matrix <- cor(data, method = "spearman")
          if (is.null(adj_matrix)) stop("Spearman returned NULL")
          networks$spearman <- construct_network("spearman", adj_matrix, use_quantile = TRUE)
        }, error = function(e) {
          errors$spearman <- e$message
          writeLines(paste("Error in Spearman:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing Spearman correlation network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (bicor$enabled) {
        start_time <- Sys.time()
        tryCatch({
          adj_matrix <- WGCNA::bicor(data)
          if (is.null(adj_matrix)) stop("Bicor returned NULL")
          networks$bicor <- construct_network("bicor", adj_matrix, use_quantile = TRUE)
        }, error = function(e) {
          errors$bicor <- e$message
          writeLines(paste("Error in Bicor:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing Biweight Midcorrelation network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (sparcc$enabled) {
        start_time <- Sys.time()
        tryCatch({
          if (!quantitative) {
            result <- do.call(SparCC.frac, c(list(data), sparcc$params[-which(names(sparcc$params) == "imax")]))
          } else {
            result <- do.call(SparCC.count, c(list(data), sparcc$params))
          }
          adj_matrix <- result$cor.w
          if (is.null(adj_matrix)) stop("SparCC returned NULL")
          networks$sparcc <- construct_network("sparcc", adj_matrix, use_quantile = TRUE)
        }, error = function(e) {
          errors$sparcc <- e$message
          writeLines(paste("Error in SparCC:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing SparCC network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (spiecEasi_mb$enabled) {
        start_time <- Sys.time()
        tryCatch({
          result <- do.call(spiec.easi, c(list(data), spiecEasi_mb$params))
          beta_matrix <- getOptBeta(result)
          inverse_cov_matrix <- SpiecEasi::symBeta(beta_matrix, mode='maxabs')
          adj_matrix <- as.matrix(inverse_cov_matrix)
          binary_matrix <- as.matrix(result$refit$stars)
          # Save network and binary network for SpiecEasi (mb)
          write.csv(adj_matrix, file.path("Network", "spiecEasi_mb.csv"), row.names = TRUE)
          write.csv(binary_matrix, file.path("Binary_Network", "spiecEasi_mb_binary.csv"), row.names = TRUE)
          networks$spiecEasi_mb <- list(adj_matrix = adj_matrix, binary_matrix = binary_matrix)
        }, error = function(e) {
          errors$spiecEasi_mb <- e$message
          writeLines(paste("Error in SpiecEasi (mb):", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing SpiecEasi (mb) network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (spiecEasi_glasso$enabled) {
        start_time <- Sys.time()
        tryCatch({
          result <- do.call(spiec.easi, c(list(data), spiecEasi_glasso$params))
          inverse_cov_matrix <- SpiecEasi::getOptCov(result)
          adj_matrix <- as.matrix(inverse_cov_matrix)
          binary_matrix <- as.matrix(result$refit$stars)
          # Save network and binary network for SpiecEasi (glasso)
          write.csv(adj_matrix, file.path("Network", "spiecEasi_glasso.csv"), row.names = TRUE)
          write.csv(binary_matrix, file.path("Binary_Network", "spiecEasi_glasso_binary.csv"), row.names = TRUE)
          networks$spiecEasi_glasso <- list(adj_matrix = adj_matrix, binary_matrix = binary_matrix)
        }, error = function(e) {
          errors$spiecEasi_glasso <- e$message
          writeLines(paste("Error in SpiecEasi (glasso):", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing SpiecEasi (glasso) network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (spring$enabled) {
        start_time <- Sys.time()
        tryCatch({
          spring_params <- spring$params
          if (!quantitative) {
            spring_params$quantitative <- FALSE
          }
          result <- do.call(SPRING, c(list(data), spring_params))
          opt.K <- result$output$stars$opt.index
          adj.K <- as.matrix(result$fit$est$path[[opt.K]])
          pcor.K <- as.matrix(SpiecEasi::symBeta(result$output$est$beta[[opt.K]], mode = 'maxabs'))
          # Save network and binary network for SPRING
          write.csv(pcor.K, file.path("Network", "spring.csv"), row.names = TRUE)
          write.csv(adj.K, file.path("Binary_Network", "spring_binary.csv"), row.names = TRUE)
          networks$spring <- list(adj_matrix = pcor.K, binary_matrix = adj.K)
        }, error = function(e) {
          errors$spring <- e$message
          writeLines(paste("Error in SPRING:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing SPRING network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (gcoda$enabled) {
        start_time <- Sys.time()
        tryCatch({
          gcoda_params <- gcoda$params
          if (quantitative) {
            gcoda_params$counts <- TRUE
          } else {
            gcoda_params$counts <- FALSE
          }
          gcoda_fun <- match.fun("gcoda")
          result <- do.call(gcoda_fun, c(list(x = data), gcoda_params))
          if (is.null(result)) stop("GCODA returned NULL")
          adj_matrix <- as.matrix(result$opt.icov)
          binary_matrix <- as.matrix(result$refit)
          # Save network and binary network for GCODA
          write.csv(adj_matrix, file.path("Network", "gcoda.csv"), row.names = TRUE)
          write.csv(binary_matrix, file.path("Binary_Network", "gcoda_binary.csv"), row.names = TRUE)
          networks$gcoda <- list(adj_matrix = adj_matrix, binary_matrix = binary_matrix)
        }, error = function(e) {
          errors$gcoda <- e$message
          writeLines(paste("Error in GCODA:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing GCODA network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (c_MI$enabled) {
        start_time <- Sys.time()
        tryCatch({
          c_MI_params <- c_MI$params
          if (quantitative) {
            c_MI_params$quantitative <- TRUE
          } else {
            c_MI_params$quantitative <- FALSE
          }
          c_MI_fun <- match.fun("conditional_MI")
          result <- do.call(c_MI_fun, c(list(data), c_MI_params))
          if (is.null(result)) stop("c_MI returned NULL")
          # Handle NaNs in Gval before proceeding
          result$Gval_order1[is.na(result$Gval_order1)] <- 0
          adj_matrix <- as.matrix(result$G_order1)
          networks$c_MI <- construct_network("c_MI", adj_matrix, use_quantile = FALSE)
        }, error = function(e) {
          errors$c_MI <- e$message
          writeLines(paste("Error in c_MI:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing c_MI network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }

      if (cclasso$enabled) {
        start_time <- Sys.time()
        tryCatch({
          cclasso_params <- cclasso$params
          if (quantitative) {
            cclasso_params$counts <- TRUE
          } else {
            cclasso_params$counts <- FALSE
          }
          cclasso_fun <- match.fun("cclasso")
          result <- do.call(cclasso_fun, c(list(x = data), cclasso_params))
          if (is.null(result)) stop("CCLasso returned NULL")
          adj_matrix <- as.matrix(result$cor_w)
          adj_matrix[is.na(adj_matrix)] <- 0
          networks$cclasso <- construct_network("cclasso", adj_matrix, use_quantile = TRUE)
        }, error = function(e) {
          errors$cclasso <- e$message
          writeLines(paste("Error in CCLasso:", e$message), con = error_log)
        })
        end_time <- Sys.time()
        cat("Constructing CCLasso network completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")
      }
    }
  ))

  # Combine binary matrices to create a final weighted network
  weighted_matrix <- matrix(0, nrow = num_features, ncol = num_features)
  rownames(weighted_matrix) <- colnames(data)
  colnames(weighted_matrix) <- colnames(data)
  methods_confirmed <- matrix(vector("list", num_features * num_features), nrow = num_features, ncol = num_features)
  rownames(methods_confirmed) <- colnames(data)
  colnames(methods_confirmed) <- colnames(data)

  for (method in names(networks)) {
    binary_matrix <- networks[[method]]$binary_matrix
    if (!is.null(binary_matrix) && nrow(binary_matrix) == num_features && ncol(binary_matrix) == num_features) {
      for (i in 1:num_features) {
        for (j in 1:num_features) {
          if (!is.na(binary_matrix[i, j]) && binary_matrix[i, j] == 1) {
            weighted_matrix[i, j] <- weighted_matrix[i, j] + 1
            methods_confirmed[[i, j]] <- c(methods_confirmed[[i, j]], method)
          }
        }
      }
    }
  }

  # Save the weighted network matrix as CSV
  diag(weighted_matrix) <- 0
  write.csv(weighted_matrix, file.path("Network", "weighted_network.csv"), row.names = TRUE)

  # Create edge lists
  edge_list <- data.frame()
  for (i in 1:num_features) {
    for (j in 1:num_features) {
      if (weighted_matrix[i, j] > 0) {
        edge_list <- rbind(edge_list, data.frame(From = rownames(weighted_matrix)[i],
                                                 To = colnames(weighted_matrix)[j],
                                                 Weight = weighted_matrix[i, j],
                                                 Methods = paste(unlist(methods_confirmed[[i, j]]), collapse = ", ")))
      }
    }
  }


  # Close the error log file
  close(error_log)

  return(list(weighted_network = weighted_matrix,
              edge_list= edge_list,
              errors = errors))
}
