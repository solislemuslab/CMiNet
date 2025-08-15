# Contributing to CMiNet

The following guidelines are designed for contributors to **CMiNet**.

## Reporting Issues and Questions

For reporting a bug, a failed function, or requesting a new feature, please open an issue in the [GitHub issue tracker](https://github.com/solislemuslab/CMiNet/issues).  
First, search through existing issues (open or closed) that might already address your question.  

If you are reporting a bug, please include:
- A minimal reproducible example (code + sample data if possible)
- The expected behavior
- The actual behavior observed
- Your R version and operating system
- Any error messages

For general questions, please:
- Check the [User Guide and Examples](https://github.com/solislemuslab/CMiNet)
- Review the online Shiny app at [https://cminet.wid.wisc.edu](https://cminet.wid.wisc.edu)  
If you still cannot find an answer, open a new GitHub issue. We will respond as soon as possible, though our team is small and we appreciate your patience.

---

## Contributing Code

To make contributions to CMiNet, you need to have a GitHub account and submit your change(s) via a **pull request** against the `development` branch from a non-`development` branch in your fork.  
Using a separate branch will make it easier to collaborate and review changes.

## Add a new algorithm to **CMiNet**.
Add a new algorithm by creating a small wrapper in R/ (e.g., R/run_newalgo.R) modeled on existing files like gcoda.R or cclasso.R; the function should accept CMiNet’s standard inputs (data, quantitative, ...) and return a symmetric taxa×taxa adjacency matrix with a zero diagonal and taxa names. Then register it in R/CMiNet.R by (i) adding newalgo = list(enabled = FALSE, params = list()) to the CMiNet() argument list and (ii) calling run_newalgo(...) inside CMiNet() and storing the result as method_results[["newalgo"]] (plus a short metadata entry for category/citation). If your method needs a package, add it to DESCRIPTION under Imports; optionally expose parameters in the Shiny app (follow the patterns in Vis_FinalNet.R / visualization.R). Please include roxygen2 docs with @references, a minimal tests/testthat/test-newalgo.R (matrix shape, symmetry, zero diagonal), run devtools::document() and devtools::check(), and submit a pull request.

### Steps for contributing:
1. **Fork** the CMiNet repository to your GitHub account.  
2. **Clone** your fork locally:  
   ```bash
   git clone https://github.com/YOUR-USERNAME/CMiNet.git
