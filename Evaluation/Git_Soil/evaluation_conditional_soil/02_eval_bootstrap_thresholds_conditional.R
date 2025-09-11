#rm(list=ls())
set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)
source("required_function.R")
count_edges <- function(network) {
  return(sum(network))
}
# Load the saved weighted networks
load("weighted_networks_list_conditional.RData")

# -------------------- LOAD & FILTER DATA --------------------
# Load your data
data <- read.csv("data/da_class_count_scabpit.csv")
dim(data)
# Filter samples
nf <- 20  # minimum number of non-zero taxa & sum threshold
taxa_data <- data[, -ncol(data)]

# Filter samples with >= nf non-zero taxa & total abundance >= nf
nonzero_taxa_counts <- rowSums(taxa_data > 0)
keep_rows <- nonzero_taxa_counts >= nf & rowSums(taxa_data) >= nf
taxa_data_filtered <- taxa_data[keep_rows, ]
dim(taxa_data_filtered)
y = taxa_data_filtered
# ------------------------------------------------------------
n_bootstrap = 50

sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
spiecEasi_mb_params= list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
spiecEasi_glasso_params =params = list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
c_MI_params = list(quantitative = TRUE, q1 = 0.7, q2 = 0.95)
cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)



# Function to compute binary network based on a given threshold
compute_binary_networks <- function(score) {
  Result_final <- data.frame(
    number_edges_cminet = numeric(),
    number_edges_phy = numeric(),
    TP = numeric(),
    F_score = numeric()
  )
  
  result <- CMiNet(
    y,
    quantitative = TRUE,
    TT = 0.95,
    pearson = list(enabled = FALSE),
    spearman = list(enabled = FALSE),
    bicor = list(enabled = FALSE),
    sparcc = list(enabled = FALSE,params=sparcc_params),
    spiecEasi_mb = list(enabled = TRUE,params = spiecEasi_mb_params),
    spiecEasi_glasso = list(enabled = TRUE,params = spiecEasi_glasso_params),
    spring = list(enabled = TRUE,params = spring_params),
    gcoda = list(enabled =TRUE, params =gcoda_params),
    c_MI  = list(enabled =TRUE,params=c_MI_params),
    cclasso = list(enabled = FALSE,params=cclasso_params)
  )
  WN = result$weighted_network
  network_final <- ifelse(WN > score, 1, 0)
  diag(network_final) = 0
  network_final[lower.tri(network_final)] <- 0  # Only upper triangle will have non-zero values
 
  
  
  for (i in 1:n_bootstrap) {
    binary_matrix <- ifelse(weighted_networks_list[[i]] > score, 1, 0)
    binary_matrix[lower.tri(binary_matrix)] <- 0  # Upper triangle only
    diag(binary_matrix) = 0
    
    binary_edges = count_edges(binary_matrix)
    # Evaluation
    EV = evaluate(dim(binary_matrix)[2], binary_matrix, network_final)
    TP <- EV[[1]][1, 1]
    F_score <- EV[[2]][1, 5]
    
    # Save results
    Result_final <- rbind(Result_final, data.frame(
      number_edges_cminet = count_edges(network_final),
      number_edges_phy = binary_edges,
      TP = TP,
      F_score = F_score
    ))
  }
  
  # Sort by F-score
  #Result_final <- Result_final[order(-Result_final$F_score), ]
  
  # Save results to CSV
  output_filename <- paste0("Result_soil_conditional_", score, ".csv")
  write.csv(Result_final, output_filename, row.names = FALSE)
}

for (score in c(0,1,2,3,4)) {
  compute_binary_networks(score)
}


