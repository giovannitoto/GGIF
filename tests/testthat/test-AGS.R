test_that("Adaptive Gibbs Sampler: wrong input parameters", {
  # data
  data("W.rda")
  data("X.rda")
  data("Y.rda")
  # correct input parameters
  seed <- 28
  stdx <- TRUE
  stdw <- TRUE
  WFormula <- formula("~ .")
  XFormula <- formula("~ .")
  kinit <- NULL
  kmax <- NULL
  kval <- 6
  nrun <- 5
  burn <- round(nrun / 4)
  thin <- 1
  start_adapt <- 5
  b0 <- 1
  b1 <- 5 * 10^(-4)
  sd_b <- 1
  sd_mu <- 1
  sd_beta <- 1
  a_theta <- 1
  b_theta <- 1
  as <- 1
  bs <- 1
  p_constant <- NULL
  alpha <- 5
  y_max <- Inf
  output <- ""
  verbose <- TRUE
  # wrong arguments
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = c(1, 2),
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'seed' not valid: it must be an integer.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = 28.01,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'seed' not valid: it must be an integer.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = "false"), "'verbose' not valid: it must be 'TRUE' or 'FALSE'.")
  expect_error(AGS_SIS(Y = 1:15, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'Y' not valid: it must be a matrix.")
  expect_error(AGS_SIS(Y = Y, X = matrix(1, nrow=4, ncol=10), W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'Y' and 'X' not compatible: the number of columns of 'Y' must be equal to the number of rows of 'X'.")
  expect_error(AGS_SIS(Y = Y, X = X, W = matrix(1, nrow=4, ncol=10),
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'Y' and 'W' not compatible: they must have the same number of rows.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = 25, kmax = 22, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'kmax' not valid: it must be an integer greater than or equal to 'kinit'.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = 30, burn = 30, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'burn' not valid: it must be an integer less than 'nrun'.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = 30, burn = 25, thin = 10, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'thin' not valid: it must be an integer less than or equal to 'nrun - burn'.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = -0.1, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = p_constant, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose), "'b0' not valid: it must be greater than or equal to 0.")
  expect_error(AGS_SIS(Y = Y, X = X, W = W,
                       seed = seed,
                       stdx = stdx, stdw = stdw,
                       WFormula = WFormula,
                       XFormula = XFormula,
                       kinit = kinit, kmax = kmax, kval = kval,
                       nrun = nrun, burn = burn, thin = thin, start_adapt = start_adapt,
                       b0 = b0, b1 = b1,
                       sd_b = sd_b, sd_mu = sd_mu, sd_beta = sd_beta,
                       a_theta = a_theta, b_theta = b_theta,
                       as = as, bs = bs,
                       p_constant = 0, alpha = alpha,
                       y_max = y_max, output = output,
                       verbose = verbose))
})
