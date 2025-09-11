rm(list = ls())
set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)

###Default parameter
sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
spiecEasi_mb_params= list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
spiecEasi_glasso_params =list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
c_MI_params = list(quantitative = TRUE, q1 = 0.9, q2 = 0.95)
cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)

Synth_SP <- function(data, ADJ, numberseed){
  n <- nrow(data)
  Cor1 <- Cor_AdjGraph(ADJ)
  X1_count <- SPRING::synthData_from_ecdf(data, Sigma = Cor1, n = n, seed = numberseed)
  return(X1_count)
}

#' Convert Adjacency Matrix to Correlation Matrix
Cor_AdjGraph <- function(ADJ){
  class(ADJ) <- "graph"
  Prec_ADJ  <- SpiecEasi::graph2prec(ADJ)
  Cor_ADJ <- cov2cor(SpiecEasi::prec2cov(Prec_ADJ))
  return(Cor_ADJ)
}
# ======================
# 0) Prepare AMGut-based counts and simulate
# ======================
# You already have amgut1.filt in your env; if not, load it (SpiecEasi has data)
if (!exists("amgut1.filt")) data("amgut1.filt", package = "SpiecEasi")

depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
#https://rdrr.io/github/GraceYoon/SPRING/man/synthData_from_ecdf.html
d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d
#graph <- SpiecEasi::make_graph('band', ncol(amgut1.filt), 2*ncol(amgut1.filt))
graph <- SpiecEasi::make_graph('band', ncol(amgut1.filt), 175)
X<- Synth_SP (amgut1.filt.cs, graph, 100)
colnames(X) <- paste0("OTU", seq_len(ncol(X)))
print(X)
A_true <- as.matrix(graph)
diag(A_true) <- 0 
edge_metrics <- function(A_est, A_true) {
  ut <- upper.tri(A_true)
  tp <- sum(A_est[ut] == 1 & A_true[ut] == 1)
  fp <- sum(A_est[ut] == 1 & A_true[ut] == 0)
  fn <- sum(A_est[ut] == 0 & A_true[ut] == 1)
  precision <- ifelse(tp + fp == 0, NA, tp/(tp+fp))
  recall    <- ifelse(tp + fn == 0, NA, tp/(tp+fn))
  f1        <- ifelse(is.na(precision) | is.na(recall) | (precision+recall)==0,
                      NA, 2*precision*recall/(precision+recall))
  jaccard   <- ifelse((tp + fp + fn) == 0, NA_real_, tp/(tp+fp+fn))
  data.frame(tp=tp, fp=fp, fn=fn,
             precision=precision, recall=recall, F1=f1, Jaccard=jaccard,
             check.names=FALSE)
}

data = X

  result <- CMiNet(
    data,
    quantitative = TRUE,
    TT = 0.95,
    pearson = list(enabled = TRUE),
    spearman = list(enabled = TRUE),
    bicor = list(enabled = TRUE),
    sparcc = list(enabled = TRUE, params = sparcc_params),
    spiecEasi_mb = list(enabled = TRUE, params = spiecEasi_mb_params),
    spiecEasi_glasso = list(enabled = TRUE, params = spiecEasi_glasso_params),
    spring = list(enabled = TRUE, params = spring_params),
    gcoda = list(enabled = TRUE, params = gcoda_params),
    c_MI = list(enabled = TRUE, params = c_MI_params),
    cclasso = list(enabled =TRUE, params = cclasso_params)
  )
  

  # After running CMiNet(), the adjacency matrices are written to "Binary_Network" and "Network".
  # We reload them here for evaluation against the true graph.
  if (!dir.exists("Binary_Network") | !dir.exists("Network")) {
    stop("Run CMiNet() first to generate network files.")
  }
  Net_pearson = read.csv("Binary_Network/pearson_binary.csv",row.names = 1)
  Net_spearman = read.csv("Binary_Network/spearman_binary.csv", row.names = 1)
  Net_bicor = read.csv("Binary_Network/bicor_binary.csv", row.names = 1)
  Net_sparcc = read.csv("Binary_Network/sparcc_binary.csv", row.names = 1)
  Net_cclasso = read.csv("Binary_Network/cclasso_binary.csv", row.names = 1)
  Net_spiecEasi_mb = read.csv("Binary_Network/spiecEasi_mb_binary.csv", row.names = 1)
  Net_spiecEasi_glasso = read.csv("Binary_Network/spiecEasi_glasso_binary.csv", row.names = 1)
  Net_spring = read.csv("Binary_Network/spring_binary.csv", row.names = 1)
  Net_gcoda= read.csv("Binary_Network/gcoda_binary.csv", row.names = 1)
  Net_c_MI = read.csv("Binary_Network/c_MI_binary.csv", row.names = 1)
  Net_cminet = read.csv("Network/weighted_network.csv", row.names = 1)
  diag(Net_pearson) <- 0
  diag(Net_spearman) <- 0
  diag(Net_bicor) <- 0
  diag(Net_sparcc) <- 0
  diag(Net_cclasso) <- 0
  diag(Net_spiecEasi_mb) <- 0
  diag(Net_spiecEasi_glasso) <- 0
  diag(Net_spring) <- 0
  diag(Net_gcoda) <- 0
  diag(Net_c_MI) <- 0
  diag(Net_cminet) <-0
  
  score = 8
  Net_cminet =   ifelse(Net_cminet > score, 1, 0)

write.csv(A_true,    "adj_TRUE.csv", row.names = FALSE)

metrics <- rbind(
  cbind(method = "pearson",    edge_metrics(Net_pearson,  A_true)),
  cbind(method = "spearman",       edge_metrics(Net_spearman, A_true)),
  cbind(method = "bicor",  edge_metrics(Net_bicor, A_true)),
  cbind(method = "sparcc",     edge_metrics(Net_sparcc,     A_true)),
  cbind(method = "cclasso",    edge_metrics(Net_cclasso,  A_true)),
  cbind(method = "spiecEasi_mb",       edge_metrics(Net_spiecEasi_mb, A_true)),
  cbind(method = "spiecEasi_glasso",  edge_metrics(Net_spiecEasi_glasso , A_true)),
  cbind(method = "spring",     edge_metrics(Net_spring,     A_true)),
  cbind(method = "gcoda",  edge_metrics( Net_gcoda , A_true)),
  cbind(method = "CMIMN",     edge_metrics(Net_c_MI,     A_true)),
  cbind(method = "Cminet",     edge_metrics(Net_cminet,     A_true)) 
)
write.csv(metrics, "metrics.csv", row.names = FALSE)
print(metrics)

