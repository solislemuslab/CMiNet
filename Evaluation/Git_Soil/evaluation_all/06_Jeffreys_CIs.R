# ============================================
# CMiNet bootstrap (from saved matrices only)
# - Loads original W (method counts) and bootstrap Wb_list
# - Computes C_ij, Jeffreys CIs, and applies joint rule
# - No bootstraps re-run
# ============================================

rm(list = ls()); set.seed(123)

# ----------------- USER INPUTS -----------------
m0     <- 8      # within-replicate support threshold
m_star <- 9      # frequency threshold on original data
c_star <- 0.95   # bootstrap confidence threshold

# paths to your saved objects
W0_csv_path   <- "weighted_network_all.csv"         # original weighted network (method counts)
Wb_rdata_path <- "weighted_networks_list_all.RData" # list of bootstrap weighted networks



# -------------- Helpers --------------
upper_idx <- function(p) which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)

to_edge_df <- function(M) {
  p <- nrow(M); idx <- upper_idx(p)
  if (nrow(idx) == 0) return(data.frame())
  data.frame(
    taxon1 = colnames(M)[idx[,2]],
    taxon2 = rownames(M)[idx[,1]],
    value  = M[idx],
    row.names = NULL, check.names = FALSE
  )
}

jeffreys_ci <- function(k, B, alpha = 0.05) {
  c(qbeta(alpha/2, k + 0.5, B - k + 0.5),
    qbeta(1 - alpha/2, k + 0.5, B - k + 0.5))
}

# -------------- 1) Load original W0 --------------
message("Loading original weighted network (method counts) from CSV...")
W0 <- as.matrix(read.csv(W0_csv_path, row.names = 1, check.names = FALSE))
storage.mode(W0) <- "numeric"
diag(W0) <- 0
W0[lower.tri(W0)] <- 0

p <- ncol(W0)
idx <- upper_idx(p)
M   <- W0[idx]  # vectorized upper triangle (original method counts)

# -------------- 2) Load bootstrap Wb_list --------------
message("Loading bootstrap weighted networks list from .RData...")
loaded_names <- load(Wb_rdata_path)   # loads objects into env; returns their names
# Try to find the bootstrap list object
if ("weighted_networks_list" %in% loaded_names) {
  Wb_list <- get("weighted_networks_list")
} else {
  # If stored under a different name, try the first list-like object
  cand_names <- loaded_names[sapply(loaded_names, function(nm) is.list(get(nm)))]
  stopifnot(length(cand_names) >= 1)
  Wb_list <- get(cand_names[1])
  message(sprintf("Using bootstrap list object: %s", cand_names[1]))
}

# Basic checks & clean-up
stopifnot(is.list(Wb_list), length(Wb_list) > 0)
B <- length(Wb_list)
for (b in seq_len(B)) {
  stopifnot(is.matrix(Wb_list[[b]]), all(dim(Wb_list[[b]]) == c(p, p)))
  diag(Wb_list[[b]]) <- 0
  Wb_list[[b]][lower.tri(Wb_list[[b]])] <- 0
}

# -------------- 3) Compute bootstrap confidence C_ij --------------
message("Computing bootstrap confidence and Jeffreys intervals from saved matrices...")

succ_vec <- numeric(length(M))   # successes across B (present given m0)
for (b in seq_len(B)) {
  mb <- Wb_list[[b]][idx]
  succ_vec <- succ_vec + as.numeric(mb >= m0)
}
C <- succ_vec / B

# Jeffreys intervals (optional, reported)
L95 <- U95 <- numeric(length(C))
for (i in seq_along(C)) {
  ci <- jeffreys_ci(succ_vec[i], B, 0.05)
  L95[i] <- ci[1]; U95[i] <- ci[2]
}

# -------------- 4) Joint rule --------------
keep_joint <- (M >= m_star) & (C >= c_star)

# -------------- 5) Export results --------------
edges <- to_edge_df(W0)  # has taxon1, taxon2, value (rename to M_ij)
colnames(edges)[3] <- "M_ij"
edges$C_ij     <- C
edges$L95_ij   <- L95
edges$U95_ij   <- U95
edges$keep_joint <- keep_joint

edges_out <- edges[order(-edges$keep_joint, -edges$M_ij, -edges$C_ij), ]

fn <- sprintf("edges_summary_from_saved_B%d_m0_%d_mstar_%d_cstar_%.2f.csv", B, m0, m_star, c_star)
write.csv(edges_out, file.path(fn), row.names = FALSE)

message("Done. Wrote: ", file.path( fn))
