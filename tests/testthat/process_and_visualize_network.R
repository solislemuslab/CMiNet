# test-vis.R (in tests/testthat/)
test_that("process_and_visualize_network runs without error and handles thresholds", {
  set.seed(123)
  mock_matrix <- matrix(runif(100, 0, 5), nrow = 10, ncol = 10)
  mock_matrix[lower.tri(mock_matrix)] <- t(mock_matrix)[lower.tri(mock_matrix)]  # Make symmetric
  diag(mock_matrix) <- 0

  taxa_names <- paste0("Taxa", 1:10)
  rownames(mock_matrix) <- colnames(mock_matrix) <- taxa_names

  thresholds <- c(max(mock_matrix)-1, max(mock_matrix)-2, max(mock_matrix)-3, max(mock_matrix)-4)
  show_labels <- c(FALSE, FALSE, FALSE, FALSE)
  node_colors <- c("white", "lightyellow", "lightgreen", "lightblue")
  edge_colors <- c("blue", "#9491D9", "#332288", "purple")

  expect_silent(
    process_and_visualize_network(mock_matrix, taxa_names, thresholds, show_labels, node_colors, edge_colors)
  )
})

