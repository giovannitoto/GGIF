test_that("AGS_SIS_RC vs AGS_SIS", {
  # data
  load("C:/Users/Giovanni/Desktop/GGIF/data/W.rda")
  load("C:/Users/Giovanni/Desktop/GGIF/data/X.rda")
  load("C:/Users/Giovanni/Desktop/GGIF/data/Y.rda")
  # correct input parameters
  seed = 292
  stdx = TRUE; stdw = TRUE
  WFormula = formula("~ ."); XFormula = formula("~ .")
  kinit = NULL; kmax = NULL; kval = 6
  nrun = 5; burn = 0; thin = 1; start_adapt = 0
  b0 = 1; b1 = 5 * 10^(-4)
  sd_b = 1; sd_mu = 1; sd_beta = 1
  a_theta = 1; b_theta = 1
  as = 1; bs = 1;
  p_constant = NULL
  alpha = 5; y_max = Inf
  output = "all"
  verbose = TRUE
  # results
  out_MCMC_RC <- AGS_SIS_RC(Y = Y, X = X, W = W,
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
  out_MCMC_C <- AGS_SIS(Y = Y, X = X, W = W,
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
  # tests
  expect_equal(out_MCMC_RC$numFactors, out_MCMC_C$numFactors, ignore_attr = TRUE)
  valid_output <- c("mu", "bmu", "beta", "eta", "lambda", "sigmacol")
  for (par in valid_output) {
    for (it in 1:length(out_MCMC_RC$numFactors)) {
      expect_equal(out_MCMC_RC[[par]][[it]], out_MCMC_C[[par]][[it]], ignore_attr = TRUE)
    }
  }
})
