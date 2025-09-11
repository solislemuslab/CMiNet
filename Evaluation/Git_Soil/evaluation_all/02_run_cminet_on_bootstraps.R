rm(list=ls())
set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)
source("required_function.R")
n_bootstrap = 50
load("50_bootstrap_soil.RData")
sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
spiecEasi_mb_params= list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
spiecEasi_glasso_params =params = list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
c_MI_params = list(quantitative = TRUE, q1 = 0.7, q2 = 0.95)
cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)
# Create a list to store weighted networks
weighted_networks_list <- list()

for (i in 1:n_bootstrap) {
  bootstrap_y <- y_bootstrap[[i]]
  result <- CMiNet(
    bootstrap_y,
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
    cclasso = list(enabled = TRUE, params = cclasso_params)
  )
  
  weighted_networks_list[[i]] <- result$weighted_network
}

# Save the weighted networks for later use
save(weighted_networks_list, file = "weighted_networks_list_all.RData")




