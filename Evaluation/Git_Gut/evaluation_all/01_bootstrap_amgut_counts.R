set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)
source("required_function.R")
data = amgut1.filt
taxa_name <- matrix(0, nrow = dim(data)[2], ncol = 2)
taxa_name[, 1] <- colnames(data)        
taxa_name[, 2] <- 1:dim(data)[2]       
colnames(taxa_name) <- c("original", "figures_name")
n_bootstrap = 50
y = data
y_bootstrap <- generate_bootstrap_data(y, n_bootstrap)

save(y_bootstrap, file = "50_bootstrap_data.RData")