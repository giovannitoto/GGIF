test_that("lposterior", {
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
  output <- "all"
  verbose <- TRUE
  # AGS
  out_MCMC <- AGS_SIS(Y = Y, X = X, W = W,
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
                      verbose = verbose)
  # wrong arguments
  expect_error(lposterior(out_MCMC, frac_sampled = 0, columns = "k"),
               "'frac_sampled' not valid: it must be a number in \\(0,1].")
  expect_error(lposterior(out_MCMC, frac_sampled = 1.1, columns = "k"),
               "'frac_sampled' not valid: it must be a number in \\(0,1].")
  expect_error(lposterior(out_MCMC, samples = c(1, 2, 3, 4.1), columns = "k"),
               "'samples' not valid: it must be a vector of integers.")
  expect_error(lposterior(out_MCMC, samples = c(1, 2, 3, 1000), columns = "k"),
               "'samples' not valid: it must be a vector of integers between 1 and the number of iterations.")
  expect_error(lposterior(out_MCMC, frac_sampled = 1, columns = "kk"),
               "'columns' not valid: it must be 'k' or 'kstar'.")
  out_MCMC$mu <- NULL
  expect_error(lposterior(out_MCMC, frac_sampled = 1, columns = "k"),
               "If 'W' is stored in the MCMC output, then also 'mu', 'bmu' must be stored in the MCMC output.")
  out_MCMC$W <- NULL
  out_MCMC$hyperparameters$alpha <- NULL
  expect_error(lposterior(out_MCMC, frac_sampled = 1, columns = "k"),
               "alpha, a_theta, b_theta, sd_b, sd_mu, sd_beta, as, bs, p_constant, y_max must be stored in 'out_MCMC\\$hyperparameters'.")
})
