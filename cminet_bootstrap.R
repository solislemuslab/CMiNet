# CMiNet bootstrap with saving bootstrap datasets and CMiNet results (Wb_list)
# Minimal bootstrap workflow for CMiNet
# 1) create B bootstrap resamples and save
# 2) run CMiNet on original + bootstraps and save
# 3) compute edge reproducibility and filter by (m*, c*)
# Keeps edges iff: (M_ij >= m_star) AND (C_ij >= c_star).
# L95_ij is computed optionally (reported but NOT used for filtering).
set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)

# ----------------- USER INPUTS -----------------
X <- amgut1.filt   # replace with your data (samples x taxa)
B       <- 50     # number of bootstrap replicates
m0      <- 8      # within-replicate support threshold
m_star  <- 9      # frequency threshold on original data
c_star  <- 0.95   # bootstrap confidence threshold

out_dir <- "bootstrap_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# -----------------------------------------------

# Utilities
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

run_cminet_get_W <- function(Xmat) {
  fit <- CMiNet::CMiNet(
    Xmat,
    quantitative = TRUE,
    TT = 0.95,
    pearson = list(enabled = TRUE),
    spearman = list(enabled = TRUE),
    bicor = list(enabled = TRUE),
    sparcc = list(enabled = TRUE),
    spiecEasi_mb = list(enabled = TRUE),
    spiecEasi_glasso = list(enabled = TRUE),
    spring = list(enabled = TRUE),
    gcoda = list(enabled = TRUE),
    c_MI  = list(enabled = TRUE),
    cclasso = list(enabled = TRUE)
  )
  W <- fit$weighted_network   # method counts M_ij
  diag(W) <- 0
  W[lower.tri(W)] <- 0
  W
}

# -------- 1) Original fit: M_ij --------
message("Running CMiNet on original data...")
W0 <- run_cminet_get_W(X)
saveRDS(W0, file.path(out_dir, "W0_weighted_network.rds"))

p   <- ncol(X)
idx <- upper_idx(p)
M   <- W0[idx]                           # vectorized upper triangle

# -------- 2) Bootstrap resampling --------
message("Running ", B, " bootstrap resamples...")

n        <- nrow(X)
succ_vec <- numeric(length(M))   # number of times edge "present"
Wb_list  <- vector("list", B)    # save each bootstrap weighted network
boot_data <- vector("list", B)   # save bootstrap abundance datasets

for (b in seq_len(B)) {
  rows <- sample.int(n, replace = TRUE)
  Xb   <- X[rows, , drop = FALSE]
  Wb   <- run_cminet_get_W(Xb)
  
  mb <- Wb[idx]
  succ_vec <- succ_vec + as.numeric(mb >= m0)
  
  Wb_list[[b]]   <- Wb
  boot_data[[b]] <- Xb
}

saveRDS(Wb_list,   file.path(out_dir, "Wb_list.rds"))
saveRDS(boot_data, file.path(out_dir, "bootstrap_datasets.rds"))

# -------- 3) Bootstrap confidence C_ij --------
C <- succ_vec / B

# Optional: Jeffreys intervals
L95 <- U95 <- numeric(length(C))
for (i in seq_along(C)) {
  ci <- jeffreys_ci(succ_vec[i], B, 0.05)
  L95[i] <- ci[1]; U95[i] <- ci[2]
}

# -------- 4) Joint rule --------
keep_joint <- (M >= m_star) & (C >= c_star)

# -------- 5) Export results --------
edges <- to_edge_df(W0)                 
colnames(edges)[3] <- "M_ij"
edges$C_ij <- C
edges$L95_ij <- L95
edges$U95_ij <- U95
edges$keep_joint <- keep_joint

edges_out <- edges[order(-edges$keep_joint, -edges$M_ij, -edges$C_ij), ]
fn <- sprintf("edges_summary_B%d_m0_%d_mstar_%d_cstar_%.2f.csv", B, m0, m_star, c_star)
write.csv(edges_out, file.path(out_dir, fn), row.names = FALSE)

