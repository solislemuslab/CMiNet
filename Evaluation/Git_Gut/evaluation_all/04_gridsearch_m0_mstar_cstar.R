## Grid search over (m0_strict, m_star, c_star) for bootstrap edge confidence
## ---- Data inputs explained ----
# WN (original, non-bootstrap):
#   - Path: Network/weighted_network.csv
#   - Source: Written by CMiNet when you run it on the ORIGINAL dataset
#             (not resampled). Same CMiNet config as used elsewhere.
#   - Meaning: WN[i, j] is the number of inference methods (out of the 10
#              enabled) that support the edge {i, j} in the original data.
#              Values are integer method-counts, not probabilities.

# weighted_networks_list (bootstrap):
#   - File: weighted_networks_list_all.RData
#   - Object inside: weighted_networks_list (a list of length B)
#   - Source: Produced by your script `scripts/02_run_cminet_on_bootstraps.R`,
#             which runs CMiNet on each bootstrap dataset and saves the
#             resulting consensus method-count matrices.
#   - Meaning: For each bootstrap b, weighted_networks_list[[b]] is a p×p matrix
#              with entry [i, j] equal to the number of methods supporting
#              edge {i, j} in that bootstrap. Same dimensions as WN and same
#              method set; diagonals should be zeroed.
rm(list=ls()); set.seed(123)
library(SPRING); library(SpiecEasi); library(CMiNet)


## ----------- Grids you want to explore -----------
m0_grid   <- c(5,6,7,8,9)                 # within-bootstrap presence
m_star_g  <- 7:10                     # method-count threshold on original WN
c_star_g  <- c(0.75, 0.80, 0.85, 0.90,0.95,1)  # confidence cutoff

WN<-as.matrix(read.csv("Network/weighted_network.csv", row.names = 1, check.names = FALSE))
storage.mode(WN) <- "numeric"; diag(WN) <- 0
p <- ncol(WN)

for (ss in (1:10)){
  score = ss
  WN2 <- ifelse(WN >= score, 1, 0)
  WN2[lower.tri(WN2)] <- 0  # Upper triangle only
  diag(WN2) = 0
  print(sum(WN2))
}

# 2) Bootstraps: list of B bootstrap consensus count matrices (same dim as WN)
load("weighted_networks_list_all.RData") # -> weighted_networks_list
B <- length(weighted_networks_list)
stopifnot(B > 0, all(sapply(weighted_networks_list, function(m) all(dim(m) == c(p,p)))))

## ----------- Precompute edge indexing & bootstrap counts -----------
UT <- which(upper.tri(WN), arr.ind = TRUE)
n_edges <- nrow(UT)
edge_i <- UT[,1]; edge_j <- UT[,2]
tax_i  <- colnames(WN)[edge_i]; tax_j <- colnames(WN)[edge_j]

# Original method counts per edge
M_e <- WN[UT]  # vector length n_edges

# Build 3D array for convenience
arr_boot <- array(0, dim = c(p, p, B))
for (b in seq_len(B)) {
  Mb <- weighted_networks_list[[b]]
  diag(Mb) <- 0
  arr_boot[,,b] <- Mb
}

# Matrix of bootstrap method counts per edge: rows = edges, cols = B
Mb_mat <- t(apply(UT, 1, function(rc) arr_boot[rc[1], rc[2], ]))  # n_edges x B
Mbar   <- rowMeans(Mb_mat)

## ----------- Function to compute confidence at a given m0_strict -----------
compute_conf <- function(m0_strict) {
  Ib_mat <- (Mb_mat >= m0_strict) * 1L
  k      <- rowSums(Ib_mat)
  Ce     <- k / B
  L95    <- qbeta(0.025, k + 0.5, B - k + 0.5)
  U95    <- qbeta(0.975, k + 0.5, B - k + 0.5)
  data.frame(C_e = Ce, k = k, L95 = L95, U95 = U95)
}



## ----------- Run grid & summarize -----------
grid_res <- data.frame()
edge_tables <- list()  # optional: store selected edges per combo

for (m0 in m0_grid) {
  conf <- compute_conf(m0)  # C_e, k, L95, U95 for this m0
  
  for (ms in m_star_g) {
    n_method_only <- sum(M_e >= ms)
    
    for (cs in c_star_g) {
      keep <- (M_e >= ms) & (conf$C_e >= cs)
      n_joint <- sum(keep)
      pruned  <- n_method_only - n_joint
      
      grid_res <- rbind(grid_res, data.frame(
        m0_strict = m0, m_star = ms, c_star = cs,
        method_only = n_method_only, joint = n_joint, pruned = pruned
      ))
      
      # OPTIONAL: save selected edge table for this combo
      key <- paste0("m0_",m0,"__mstar_",ms,"__c_",cs)
      edge_tables[[key]] <- data.frame(
        taxon_i = tax_i[keep], taxon_j = tax_j[keep],
        M_e = M_e[keep], C_e = conf$C_e[keep], k = conf$k[keep],
        L95 = conf$L95[keep], U95 = conf$U95[keep], Mbar = Mbar[keep],
        stringsAsFactors = FALSE
      )
    }
  }
}

# Save grid summary
write.csv(grid_res, "grid_summary_m0_mstar_cstar.csv", row.names = FALSE)

# OPTIONAL: write a few edge lists to disk (choose combos you care about)
to_write <- c("m0_5__mstar_8__c_0.8", "m0_6__mstar_8__c_0.8", "m0_6__mstar_8__c_0.85")
for (nm in to_write) {
  if (!is.null(edge_tables[[nm]])) {
    write.csv(edge_tables[[nm]], paste0("edges_", nm, ".csv"), row.names = FALSE)
  }
}

## ----------- Handy diagnostics/plots for a chosen combo -----------
# Pick a combo to visualize (you can change these numbers quickly)
m0_vis <- 5; m_star_vis <- 8; c_star_vis <- 0.8
conf_vis <- compute_conf(m0_vis)
cand <- (M_e >= m_star_vis)

png("hist_Ce_candidates.png", 900, 700, res=150)
hist(conf_vis$C_e[cand], breaks=20,
     main=sprintf("C_e for candidate edges (M_e ≥ %d, m0=%d)", m_star_vis, m0_vis),
     xlab="C_e")
abline(v=c_star_vis, lty=2)
dev.off()

png("scatter_Me_vs_Ce.png", 900, 700, res=150)
plot(M_e, conf_vis$C_e, pch=20,
     xlab="Original Method Count (M_e)",
     ylab=sprintf("Bootstrap Confidence (C_e), m0=%d", m0_vis))
abline(v=m_star_vis, lty=2); abline(h=c_star_vis, lty=2)
dev.off()

# Print a small concise table to console
subset(grid_res, m0_strict %in% m0_grid & m_star==9 & c_star %in% c(0.8,0.85,0.9,0.95))
