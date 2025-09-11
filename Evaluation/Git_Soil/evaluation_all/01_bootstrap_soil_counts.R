set.seed(123)
library(SpiecEasi)
library(SPRING)
library(CMiNet)
source("required_function.R")
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
n_bootstrap = 50
y_bootstrap <- generate_bootstrap_data(y, n_bootstrap)

save(y_bootstrap, file = "50_bootstrap_soil.RData")