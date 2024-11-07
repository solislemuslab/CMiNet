## CMiNet: Consensus Microbiome Network Algorithm
<img src="image/logo.png" style="width:13%;" align=right>

## Description
<div align="justify">
CMiNet is an R package designed to generate consensus microbiome networks by integrating results from multiple network construction algorithms. This tool is specifically tailored for microbiome data, where capturing the intricate relationships between microbial taxa is essential to understanding complex biological systems and their impacts on health and disease.

The package employs a range of established algorithms, including Pearson and Spearman correlation, Biweight midcorrelation, Sparse Correlations for Compositional data (SparCC), Sparse InversE Covariance estimation for Ecological Association and Statistical Inference (SpiecEasi), Semi-Parametric Rank-based Correlation and Partial Correlation Estimation (SPRING), Generalized Co-Occurrence Differential Abundance analysis (gCoda), Correlation Inference for Compositional Data through Lasso (CCLasso), and a novel algorithm based on conditional mutual information (c_MI). These algorithms construct individual microbial association networks, which CMiNet then combines into a single, weighted consensus network. By leveraging the strengths of each method, CMiNet provides a comprehensive and reliable representation of microbial interactions.
</div>
<img src="image/CMiNet-Page-2.jpg" style="width:50%;" align=right>

## Table of Contents
- [Methods Included in CMiNet](#methods-included-in-cminet)
- [Installation](#installation)
- [Running CMiNet Package](#Running-CMiNet-Package)
  - [CMiNet Function](#CMiNet-Function)
  - [process_and_visualize_network Function](#process_and_visualize_network-Function)
  - [plot_hamming_distances Function](#plot_hamming_distances-Function)
  - [plot_network Function](#plot_network-Function)

## Methods Included in CMiNet
Algorithms Applied in CMiNet:
- Pearson coefficient (cor() from stats package)
- Spearman coefficient (cor() from stats package)
- Biweight Midcorrelation (bicor() from WGCNA package)
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
If there are any errors during installation, please install the missing dependencies manually.
In particular the automatic installation of SPRING and SpiecEasi (only available on GitHub) does sometimes not work. These packages can be installed as follows (the order is important because SPRING depends on SpiecEasi):
```bash
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")
```
If you have difficulty installing SpiecEasi, download the tar or zip file from this [link](https://github.com/zdk123/SpiecEasi/releases/tag/v1.1.1). Then, install it in R.

## Running CMiNet Package
<div align="justify">
The CMiNet Package contains four main functions:

1. **CMiNet**: This function constructs a consensus network from microbiome data using multiple methods.

2. **process_and_visualize_network**: This function processes a weighted microbiome network and visualizes it across different thresholds. Each threshold represents a minimum edge weight required for inclusion in the network plot.

3. **plot_hamming_distances**: Calculates the Hamming distance, common edges, and the number of edges for each pair of resulting network matrices.

4. **plot_network**: This function generates a network plot from a final network produced by CMiNet.
</div>

## CMiNet Function
### loading the Data
We use the American Gut data from [SpiecEasi package](https://github.com/zdk123/SpiecEasi) to run CMiNet algorithm to construct consensus microbiome network. 
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
We designed the package to allow users to adjust the default parameters of each algorithm according to their preferences and specific research needs.
- sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
- spiecEasi_mb_params= list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
- spiecEasi_glasso_params =params = list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
- spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
- gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
- c_MI_params = list(quantitative = TRUE, q1 = 0.7, q2 = 0.95)
- cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)
### Construct weighted network by CMiNet Function
- Description: The CMiNet function constructs a consensus network from microbiome data by integrating multiple network construction algorithms. This approach combines individual microbial association networks into a single weighted network, enhancing reliability by leveraging the strengths of each method.
- Parameters:
  - data: The input microbiome abundance matrix, where rows represent samples and columns represent taxa.
  - quantitative (default = TRUE): Set to TRUE if the data contains quantitative abundance counts; set to FALSE otherwise.
  - TT: The threshold value used for threshold-dependent algorithms (pearson, spearman, bicor, SparCC, and CCLasso), set to the 0.95 quantile by default to promote sparsity in the consensus network.
  - Algorithm-Specific Parameters:
    - pearson, spearman, bicor: Boolean flags to enable these correlation-based methods. When set to TRUE, the algorithm runs on the data; when set to FALSE, it is excluded from the analysis.
    - SparCC, SpiecEasi_mb, SpiecEasi_glasso, SPRING, gCoda, c_MI, CCLasso: Lists of parameters for each algorithm, allowing customization of settings such as number of cores (ncores), threshold values (lambda.min.ratio), and other algorithm-specific configurations.
    - Setting the enabled parameter to TRUE for an algorithm will include it in the consensus network construction, while setting it to FALSE will skip the algorithm, allowing users to select a subset of methods for analysis.
- Returns: A list containing:
  - weighted_network: A weighted network matrix representing the consensus network.
  - edge_list: A matrix with three columns—first and second columns indicate node IDs, and the third column shows the edge weight values.
  - errors: Any errors from algorithms that could not run on the data, providing feedback for troubleshooting.

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
## process_and_visualize_network Function:
- Description: The process_and_visualize_network function processes a weighted microbiome network and visualizes it across specified thresholds. Each threshold filters edges based on minimum weight, helping to analyze network structure at varying levels of connectivity.

- Parameters:
  - weighted_network: The input weighted network matrix generated by the CMiNet function.
  - taxa_names: Vector of taxa names for labeling nodes in the network.
  - thresholds: A vector of threshold values defining minimum edge weights for inclusion in the network plot.
  - show_labels: Boolean vector indicating whether to display taxa names for each threshold (e.g., TRUE or FALSE).
  - node_colors: Vector of colors for nodes in each threshold-based plot.
  - edge_colors: Vector of colors for edges in each threshold-based plot.
- Returns:
A set of network plots for each threshold, visualizing taxa relationships based on edge weights. Each plot includes details on the number of nodes and edges in the network, and also the maximum degree value .
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
<img src="image/network_9,8,7,6.jpeg" style="width:90%;">

## plot_hamming_distances Function
- Description:
The plot_hamming_distances function calculates the Hamming distance, number of common edges, and total edges for each pair of network matrices. It generates a plot to visualize these distances, providing insights into structural differences between networks generated by different algorithms or parameters.

- Parameters:
  - network_type (e.g., "Binary_Network"): The directory or folder where all the binary network matrices are saved.
  - top_n_pairs: An integer specifying the top number of network pairs to display based on the smallest Hamming distances.
  - output_filename: The filename for saving the Hamming distances plot as an image file.
- Returns:
A plot visualizing the Hamming distances between selected network pairs. The plot is optionally saved as an image file to the specified output_filename.

```bash
plot_hamming_distances("Binary_Network", top_n_pairs = 10, output_filename = "hamming_distances_plot.jpeg")
```
<img src="image/hamming_distances_plot.jpeg" style="width:90%;">

## plot_network Function
- Description:
The plot_network function generates a visual representation of the final consensus network by applying a score threshold. This threshold filters the network to display only the most significant edges, allowing users to examine the core structure of the microbiome network based on selected nodes and edges.

- Parameters:
  - network_final: A binarized network matrix derived from the weighted consensus network. The matrix is filtered based on a numeric Score threshold, where only edges with weights above this threshold are included.
  - node_color: Specifies the color for nodes in the plot.
  - edge_color: Specifies the color for edges in the plot.
  - label_color: Specifies the color for node labels in the plot.
- Returns:
A figure displaying the final microbiome network, filtered by the score threshold. Nodes, edges, and labels are color-coded according to the specified parameters, providing a clear view of significant interactions within the network.
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
<img src="image/network_plot.png" style="width:90%;">




