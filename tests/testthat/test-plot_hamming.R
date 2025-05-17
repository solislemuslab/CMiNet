test_that("plot_hamming_distances creates a plot from a folder of matrices", {
  # Create a temporary folder
  tmp_dir <- tempfile("Binary_Network")
  dir.create(tmp_dir)

  # Create mock binary matrices and save as CSV
  set.seed(123)
  for (i in 1:3) {
    mat <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
    write.csv(mat, file = file.path(tmp_dir, paste0("method", i, ".csv")), row.names = FALSE)
  }

  # Temporary output filename
  output_file <- tempfile(fileext = ".jpeg")

  # Run the function
  plot_hamming_distances(tmp_dir, top_n_pairs = 3, output_filename = output_file)

  # Check that file was created
  expect_true(file.exists(output_file))

  # Clean up
  unlink(tmp_dir, recursive = TRUE)
  file.remove(output_file)
})
