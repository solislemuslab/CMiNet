test_that("CMiNet returns a valid network result list", {
  set.seed(123)
  mock_data <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  colnames(mock_data) <- paste0("Taxa", 1:10)

  sparcc_params = list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
  spiecEasi_mb_params = list(method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 20, ncores = 4))
  spiecEasi_glasso_params = list(method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 15, pulsar.params = list(rep.num = 50))
  spring_params = list(Rmethod = "original", quantitative = TRUE, ncores = 5, lambdaseq = "data-specific", nlambda = 15, rep.num = 20)
  gcoda_params = list(counts = FALSE, pseudo = 0.5, lambda.min.ratio = 1e-4, nlambda = 15, ebic.gamma = 0.5)
  c_MI_params = list(quantitative = TRUE, q1 = 0.7, q2 = 0.95)
  cclasso_params = list(counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1), k_max = 20, n_boot = 20)

  result <- CMiNet(
    mock_data,
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
    c_MI  = list(enabled = TRUE, params = c_MI_params),
    cclasso = list(enabled = TRUE, params = cclasso_params)
  )

  expect_type(result, "list")
  expect_true(all(c("weighted_network", "edge_list", "errors") %in% names(result)))
})
