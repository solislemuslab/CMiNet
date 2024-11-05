# <span style="color:blue;">C</span>Mi<span style="color:blue;">N</span>et: <span style="color:blue;">C</span>onsensus <span style="color:blue;">M</span>icrobiome <span style="color:blue;">N</span>etwork <span style="color:blue;">A</span>lgorithm
<img src="image/logo.png" style="width:40%;" align=right>

## Description
CMiNet is an R package designed to generate consensus microbiome networks by integrating results from multiple network construction algorithms. This tool is specifically tailored for microbiome data, where capturing the intricate relationships between microbial taxa is essential to understanding complex biological systems and their impacts on health and disease.

The package employs a range of established algorithms, including Pearson and Spearman correlation, Biweight midcorrelation, Sparse Correlations for Compositional data (SparCC), Sparse InversE Covariance estimation for Ecological Association and Statistical Inference (SpiecEasi), Semi-Parametric Rank-based Correlation and Partial Correlation Estimation (SPRING), Generalized Co-Occurrence Differential Abundance analysis (GCODA), Correlation Inference for Compositional Data through Lasso (CCLasso), and a novel algorithm based on conditional mutual information. These algorithms construct individual microbial association networks, which CMiNet then combines into a single, weighted consensus network. By leveraging the strengths of each method, CMiNet provides a comprehensive and reliable representation of microbial interactions.
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
# Required packages
```bash
install.packages("devtools")
devtools::install_github("rosaaghdam/CMiNet")
```

## Running CMiNet
We use the American Gut data from [SpiecEasi package](https://github.com/zdk123/SpiecEasi)) to run CMiNet algorithm to construct consensus microbiome network. 

First, load CMiNet and the American Gut Project data (included with the [SpiecEasi package](https://github.com/zdk123/SpiecEasi)), which is automatically loaded alongside CMiNet).
```bash
library(CMiNet)
data("amgut1.filt")
```
CMiNet Package contains 4 main funcitons:
- CMiNet: This function constructs a consensus network from microbiome data using multiple methods.
- plot_hamming_distances: Calculates the Hamming distance, common edges, and number of edges
for each pair of resulted network matrices
- visualization: This function processes a weighted microbiome network and visualizes it across different thresholds. Each threshold represents a minimum edge weight required for inclusion in the network plot.
- Vis_FinalNet: This function generates a network plot from a final network resulted by CMiNet.

As an example, we run CMiNet Package on American Gut data:
  

