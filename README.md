## CMiNet: Consensus Microbiome Network Algorithm
<img src="image/logo.png" style="width:30%;" align=right>

## Description
CMiNet is an R package designed to generate consensus microbiome networks by integrating results from multiple network construction algorithms. This tool is specifically tailored for microbiome data, where capturing the intricate relationships between microbial taxa is essential to understanding complex biological systems and their impacts on health and disease.

The package employs a range of established algorithms, including Pearson and Spearman correlation, Biweight midcorrelation, Sparse Correlations for Compositional data (SparCC), Sparse InversE Covariance estimation for Ecological Association and Statistical Inference (SpiecEasi), Semi-Parametric Rank-based Correlation and Partial Correlation Estimation (SPRING), Generalized Co-Occurrence Differential Abundance analysis (GCODA), Correlation Inference for Compositional Data through Lasso (CCLasso), and a novel algorithm based on conditional mutual information. These algorithms construct individual microbial association networks, which CMiNet then combines into a single, weighted consensus network. By leveraging the strengths of each method, CMiNet provides a comprehensive and reliable representation of microbial interactions.
<img src="image/CMiNet-Page-2.jpg" style="width:50%;" align=right>

## Table of Contents
- [Methods Included in CMiNet](#methods-included-in-cminet)
- [Installation](#installation)
- [Running CMiNet](#running)

## Methods Included in CMiNet
Algorithms Applied in CMiNet:
- Pearson coefficient (cor() from stats package)
- Spearman coefficient (cor() from stats package)
- Biweight Midcorrelation bicor() from WGCNA package
- SparCC ([R code on GitHub](https://github.com/huayingfang/CCLasso/blob/master/R/SparCC.R))
- CCLasso ([R code on GitHub](https://github.com/huayingfang/CCLasso/tree/master))
- SpiecEasi ([SpiecEasi package](https://github.com/zdk123/SpiecEasi))
- SPRING ([SPRING package](https://github.com/GraceYoon/SPRING))
- CMI (Conditional Mutual Information-Based Microbiome Network Construction)
- gCoda ([R code on GitHub](https://github.com/huayingfang/gCoda))

## Installation
```bash
install.packages("devtools")
devtools::install_github("rosaaghdam/CMiNet")
```

## Running CMiNet Package

The CMiNet Package contains four main functions:

1. **CMiNet**: This function constructs a consensus network from microbiome data using multiple methods.

2. **process_and_visualize_network**: This function processes a weighted microbiome network and visualizes it across different thresholds. Each threshold represents a minimum edge weight required for inclusion in the network plot.

3. **plot_hamming_distances**: Calculates the Hamming distance, common edges, and the number of edges for each pair of resulting network matrices.

4. **Vis_FinalNet**: This function generates a network plot from a final network produced by CMiNet.

As an example, we run CMiNet Package on American Gut data:

## 1. **CMiNet** Function
### loading the Data
We use the American Gut data from [SpiecEasi package](https://github.com/zdk123/SpiecEasi)) to run CMiNet algorithm to construct consensus microbiome network. 

First, load CMiNet and the American Gut Project data (included with the [SpiecEasi package](https://github.com/zdk123/SpiecEasi)), which is automatically loaded alongside CMiNet).
```bash
library(CMiNet)
data("amgut1.filt")
taxa_name <- matrix(0, nrow = dim(data)[2], ncol = 2)
taxa_name[, 1] <- colnames(data)        
taxa_name[, 2] <- 1:dim(data)[2]       
colnames(taxa_name) <- c("original taxa", "taxa name in figures")
```
### Define the parameter on all Algorithms
We design the package that you can change the paramters of all algorithms based on your interest.
- sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
- spiecEasi_mb_params= list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
- spiecEasi_glasso_params =params = list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
- spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
- gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
- c_MI_params = list(quantitative = TRUE, q1 = 0.7, q2 = 0.95)
- cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)
### Construct weighted network by CMiNet Function
When the data is abundance matrix (Count value), set quantitative = TRUE, otherwise quantitative = FALSE.
If you are interested to use some of these algorithms change TRUE to FALSE.
TT is the value for thresold-depend algortihms, since we suppose that the microbiome network is sparce this value set as quantime 0.95. 
```bash
result <- CMiNet(
  data,
  quantitative = TRUE,
  TT = 0.95,
  pearson = list(enabled = TRUE),
  spearman = list(enabled = TRUE),
  bicor = list(enabled = TRUE),
  sparcc = list(enabled = TRUE,params=sparcc_params),
  spiecEasi_mb = list(enabled = TRUE,params = spiecEasi_mb_params),
  spiecEasi_glasso = list(enabled = TRUE,params = spiecEasi_glasso_params),
  spring = list(enabled = TRUE,params = spring_params),
  gcoda = list(enabled =TRUE, params =gcoda_params),
  c_MI  = list(enabled =TRUE,params=c_MI_params),
  cclasso = list(enabled = TRUE,params=cclasso_params)
)
```
## 2. **process_and_visualize_network** Function:
We plotted weighted_network resulted by CMiNet Function based on 4 thresolds:

```bash
weighted_network = result$weighted_network
taxa_names <- taxa_name[, 2]
rownames(weighted_network) = colnames(weighted_network)=taxa_names
thresholds <- c(max(weighted_network)-1,max(weighted_network)-2,max(weighted_network)-3,max(weighted_network)-4)  # Thresholds for visualization
show_labels <- c(FALSE, FALSE, FALSE, FALSE)  # Show taxa names (TRUE/FALSE) for each threshold
node_colors <- c("white", "lightyellow",  "lightgreen","lightblue")  # Node colors for each threshold
edge_colors <- c("blue", "#9491D9", "#332288", "purple")  # Edge colors for each threshold
process_and_visualize_network(weighted_network, taxa_names, thresholds, show_labels, node_colors, edge_colors)
```
<img src="image/network_9,8,7,6.jpeg" style="width:60%;">
## 3. **plot_hamming_distances** Function
```bash
plot_hamming_distances("Binary_Network", top_n_pairs = 10, output_filename = "hamming_distances_plot.jpeg")
```
<img src="image/hamming_distances_plot.jpeg" style="width:60%;">

## 4. **plot_network** Function
```bash
score = 7
WN = result$weighted_network
taxa_names <- taxa_name[, 2]
rownames(WN) = colnames(WN)=taxa_names
network_final <- ifelse(WN > score, 1, 0)
network_final[lower.tri(network_final)] <- 0  # Only upper triangle will have non-zero values
sum(network_final)
plot_network(network_final, node_color = "skyblue", edge_color = "grey", label_color = "black")
```
<img src="image/network_plot.png" style="width:60%;">




