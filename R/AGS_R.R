#' AGS for GIF models with SIS prior (R)
#'
#' Implementation in R of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
#'
#' @param Y A nxp matrix of counts.
#' @param X_mean A pxq matrix of meta-covariates; the variables must be numeric or factors.
#' @param X_cov A pxq matrix of meta-covariates; the variables must be numeric or factors.
#' @param W A nxc matrix of covariates; the variables must be numeric or factors.
#' @param seed Seed. Default is 28.
#' @param stdx Boolean: if TRUE, numeric meta-covariates are standardized.
#' @param stdw Boolean: if TRUE, numeric covariates are standardized.
#' @param XmeanFormula Formula specifying  the meta-covariates inserted in the model.
#' @param XcovFormula Formula specifying  the meta-covariates inserted in the model.
#' @param WFormula Formula specifying  the covariates inserted in the model.
#' @param kinit An integer minimun number of latent factors. Default is \code{min(floor(log(p)*kval),p)}.
#' @param kmax Maximum number of latent factors. Default is \code{p+1}.
#' @param kval An integer initial number of latent factors. Default is 6.
#' @param nrun An integer number of iterations. Default is 100.
#' @param burn An integer number of burning iterations (number of iterations to discard). Default is \code{round(nrun/4)}.
#' @param thin An integer thinning value (number of iterations to skip between saving iterations). Default is 1.
#' @param start_adapt An integer number of iterations before adaptation. Default is 50.
#' @param y_max A fixed and known upper bound for the values in \code{Y}. Default is \code{Inf}.
# Parameters for SIS:
#' @param b0,b1 Positive constants for the adaptive probability. Default is \code{c(1,5*10^(-4))}.
#' @param sd_b Standard deviation for \eqn{b}. Default is 1.
#' @param sd_mu Standard deviation for \eqn{\mu}. Default is 1.
#' @param sd_beta Standard deviation for \eqn{\beta}. Default is 1.
#' @param a_theta,b_theta Default is \code{c(1,1)}.
#' @param as,bs \eqn{\sigma^{-2} \sim Gamma(\code{as}, \code{bs})}. Default is \code{c(1,1)}.
#' @param p_constant Factor probability constant. Default is \code{10*exp(1)*log(p)/p}.
#' @param alpha \eqn{Beta(1, \code{alpha})} for pi.
# Other settings:
#' @param output Default is \code{"all"}.
#' @param verbose If TRUE, print the number of active factors every 10 iterations. Default is \code{TRUE}.
#'
#' @return A list.
#'
#' @seealso \code{\link{update_b_mu_R}}, \code{\link{update_mu_R}}, \code{\link{update_eta_R}}, \code{\link{update_Lambda_star_R}}, \code{\link{update_d_R}}, \code{\link{update_Phi_R}}
#'
#' @importFrom stats formula model.matrix plogis rbeta rbinom rgamma rnorm runif
#'
#' @export
AGS_SIS_R <- function(Y,
                      X_mean = NULL, X_cov = NULL, W = NULL,
                      seed = 28,
                      stdx = TRUE, stdw = TRUE,
                      WFormula = formula("~ ."),
                      XmeanFormula = formula("~ ."),
                      XcovFormula = formula("~ ."),
                      kinit = NULL, kmax = NULL, kval = 6,
                      nrun = 100, burn = round(nrun/4), thin = 1,
                      start_adapt = 50,
                      b0 = 1, b1 = 5*10^(-4),
                      sd_b = 1, sd_mu = 1, sd_beta = 1,
                      a_theta = 1, b_theta = 1,
                      as = 1, bs = 1,
                      p_constant = NULL, alpha,
                      y_max = Inf,
                      output = "all",
                      verbose = TRUE) {
  # -------------------------------------------------------------------------- #
  # set seed
  if((length(seed) != 1) || !is.numeric(seed) || (seed != round(seed))) {
    stop("'seed' not valid: it must be an integer.")
  } else {
    set.seed(seed)
  }
  # -------------------------------------------------------------------------- #
  if((length(verbose) != 1) || !is.logical(verbose)) {
    stop("'verbose' not valid: it must be 'TRUE' or 'FALSE'.")
  }
  # -------------------------------------------------------------------------- #
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  # -------------------------------------------------------------------------- #
  if(is.null(X_mean)) {
    X_mean <- matrix(1, nrow = p, ncol = 1)
  } else {
    if(!is.data.frame(X_mean) && !is.matrix(X_mean)) {
      stop("'X_mean' not valid: it must be a matrix or a dataframe.")
    }
    if(p != nrow(X_mean)) {
      stop("'Y' and 'X_mean' not compatible: the number of columns of 'Y' must be equal to the number of rows of 'X_mean'.")
    }
    if((length(stdx) != 1) || !is.logical(stdx)) {
      stop("'stdx' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdx & is.data.frame(X_mean)) {
      is.fact.xm = sapply(X_mean, is.factor)
      X_mean[, is.fact.xm == FALSE] <- scale(X_mean[, is.fact.xm == FALSE])
      X_mean <- model.matrix(XmeanFormula, X_mean)
    }
  }
  if(is.null(X_cov)) {
    X_cov <- X_mean
  } else {
    if(!is.data.frame(X_cov) && !is.matrix(X_cov)) {
      stop("'X_cov' not valid: it must be a matrix or a dataframe.")
    }
    if(p != nrow(X_cov)) {
      stop("'Y' and 'X_cov' not compatible: the number of columns of 'Y' must be equal to the number of rows of 'X_cov'.")
    }
    if((length(stdx) != 1) || !is.logical(stdx)) {
      stop("'stdx' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdx & is.data.frame(X_cov)) {
      is.fact.xc <- sapply(X_cov, is.factor)
      X_cov[, is.fact.xc == FALSE] <- scale(X_cov[, is.fact.xc == FALSE])
      X_cov <- model.matrix(XcovFormula, X_cov)
    }
  }
  q_mean <- ncol(X_mean)
  q_cov <- ncol(X_cov)
  # -------------------------------------------------------------------------- #
  Wnull <- is.null(W)
  if(!Wnull) {
    if(!is.data.frame(W) && !is.matrix(W)) {
      stop("'W' not valid: it must be a matrix or a dataframe.")
    }
    if(n != nrow(W)) {
      stop("'Y' and 'W' not compatible: they must have the same number of rows.")
    }
    if((length(stdw) != 1) || !is.logical(stdw)) {
      stop("'stdw' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdw & is.data.frame(W)) {
      is.fact.w <- sapply(W, is.factor)
      W[, is.fact.w == FALSE] <-  scale(W[, is.fact.w == FALSE])
      W <- model.matrix(WFormula, W)
    }
    c <- ncol(W)
  }
  # -------------------------------------------------------------------------- #
  if((length(kval) != 1) || !is.numeric(kval) || (kval != round(kval))) {
    stop("'kval' not valid: it must be an integer.")
  }
  if(is.null(kinit)) {
    kinit <- min(floor(log(p) * kval), p)
  } else if((length(kinit) != 1) || !is.numeric(kinit) || (kinit != round(kinit))) {
    stop("'kinit' not valid: it must be an integer.")
  }
  if(is.null(kmax)) {
    kmax <- p + 1
  } else if((length(kmax) != 1) || !is.numeric(kmax) || (kmax != round(kmax)) || (kmax < kinit)) {
    stop("'kmax' not valid: it must be an integer greater than or equal to 'kinit'.")
  }
  k <- kinit       # number of factors to start with (active and inactive)
  kstar <- k - 1   # number of active factors
  # -------------------------------------------------------------------------- #
  if((length(nrun) != 1) || !is.numeric(nrun) || (nrun != round(nrun))) {
    stop("'nrun' not valid: it must be an integer.")
  }
  if((length(burn) != 1) || !is.numeric(burn) || (burn != round(burn)) || (burn >= nrun)) {
    stop("'burn' not valid: it must be an integer less than 'nrun'.")
  }
  if((length(thin) != 1) || !is.numeric(thin) || (thin != round(thin)) || (thin > nrun - burn)) {
    stop("'thin' not valid: it must be an integer less than or equal to 'nrun - burn'.")
  }
  if((length(start_adapt) != 1) || !is.numeric(start_adapt) || (start_adapt != round(start_adapt))) {
    stop("'start_adapt' not valid: it must be an integer.")
  }
  # number of posterior samples
  sp <- floor((nrun - burn) / thin)
  # -------------------------------------------------------------------------- #
  if((length(y_max) != 1) || !is.numeric(y_max)) {
    stop("'y_max' not valid: it must be a number or 'Inf'.")
  }
  # Transformation for the response
  a_j <- function(j, y_max) {
    val <- j
    val[j == y_max + 1] <- Inf
    val
  }
  # Bounds for truncated normal
  a_y <- a_yp1 <- matrix(NA, nrow = n, ncol = p)
  for(j in 1:p) {
    a_y[, j] <- a_j(Y[, j], y_max)        # a_y = a_j(y)
    a_yp1[, j] <- a_j(Y[, j] + 1, y_max)  # a_yp1 = a_j(y + 1)
  }
  # Replace NA with 0/Inf in a_y/a_yp1
  a_y[is.na(Y)] <- 0      # log(0)=-Inf
  a_yp1[is.na(Y)] <- Inf  # log(Inf)=Inf
  # -------------------------------------------------------------------------- #
  # Adaptive probability
  if((length(b0) != 1) || !is.numeric(b0) || (b0 < 0)) {
    stop("'b0' not valid: it must be greater than or equal to 0.")
  }
  if((length(b1) != 1) || !is.numeric(b1) || (b1 <= 0)) {
    stop("'b1' not valid: it must be greater than 0.")
  }
  prob <- 1 / exp(b0 + b1 * seq(1, nrun))
  uu <- runif(nrun)
  # -------------------------------------------------------------------------- #
  # p constant (c_p)
  if(is.null(p_constant)) {
    p_constant <- 10 * exp(1) * log(p) / p
  } else if((length(p_constant) != 1) || !is.numeric(p_constant) || (p_constant <= 0) || (p_constant >= 1)) {
    stop("'p_constant' not valid: it must be a number in (0,1).")
  }
  # -------------------------------------------------------------------------- #
  # check fixed parameters
  if((length(sd_b) != 1) || !is.numeric(sd_b) || (sd_b <= 0)) {
    stop("'sd_b' not valid: it must be greater than 0.")
  }
  if((length(sd_mu) != 1) || !is.numeric(sd_mu) || (sd_mu <= 0)) {
    stop("'sd_mu' not valid: it must be greater than 0.")
  }
  if((length(sd_beta) != 1) || !is.numeric(sd_beta) || (sd_beta <= 0)) {
    stop("'sd_beta' not valid: it must be greater than 0.")
  }
  if((length(a_theta) != 1) || !is.numeric(a_theta) || (a_theta <= 0)) {
    stop("'a_theta' not valid: it must be greater than 0.")
  }
  if((length(b_theta) != 1) || !is.numeric(b_theta) || (b_theta <= 0)) {
    stop("'b_theta' not valid: it must be greater than 0.")
  }
  if((length(as) != 1) || !is.numeric(as) || (as <= 0)) {
    stop("'as' not valid: it must be greater than 0.")
  }
  if((length(bs) != 1) || !is.numeric(bs) || (bs <= 0)) {
    stop("'bs' not valid: it must be greater than 0.")
  }
  if((length(alpha) != 1) || !is.numeric(alpha) || (alpha < 0)) {
    stop("'alpha' not valid: it must be greater than or equal to 0.")
  }
  # -------------------------------------------------------------------------- #
  # Initialize sigma^-2
  ps <- rgamma(p, as, bs)
  # Initialize the nxp matrix of latent normal variables
  # Z <- matrix(NA, nrow = n, ncol = p)
  # Initialize parameters related to the covariates, if W exists
  if(!Wnull) {
    mu <- matrix(rnorm(c * p, 0, sd_mu), nrow = p, ncol = c)    # mean coeff of the data
    b_mu <- matrix(rnorm(q_mean * c), nrow = q_mean, ncol = c)  # x effects on mu coeff
    # precision of mu and b_mu
    prec_b  <- 1 / (sd_b)  ^ 2
    prec_mu <- 1 / (sd_mu) ^ 2
  }
  # Initialize lambda star (pxq)
  Lambda_star <- matrix(rnorm(p * k), nrow = p, ncol = k)   # loading matrix
  # Initialize eta (nxk)
  eta <- matrix(rnorm(n * k), nrow = n, ncol = k)           # latent factors
  # Initialize Beta (qxk) and pred (pxk)
  Beta <- matrix(rnorm(q_cov * k), nrow = q_cov, ncol = k)  # traits effect on local shrinkage
  pred <- X_cov %*% Beta                                    # local shrinkage coefficients
  logit <- plogis(pred)
  # Initialize Phi pxk
  Phi <- matrix(rbinom(p * k, size = 1, prob = p_constant), nrow = p, ncol = k)
  # Initialize pi_h, h = 1, ..., k
  v <- c(rbeta(k - 1, shape1 = 1, shape2 = alpha), 1)
  w <- v * c(1, cumprod(1 - v[-k]))  # product up to  l - 1
  d <- rep(k, k)                 # augmented data
  rho <- rep(1, k)                   # preallocation for Bernoulli
  # Initialize the precision matrix of lambda star
  Plam <- diag(rgamma(k, a_theta, b_theta))
  # Compute Lambda (pxk)
  Lambda <- t(t(Lambda_star) * sqrt(rho)) * sqrt(Phi)
  # -------------------------------------------------------------------------- #
  # Allocate output object memory
  valid_outputs <- c("mu",          # coefSamples          : pxc
                     "bmu",         # bmuCoefSamples       : qxc
                     "beta",        # shrinkCoefSamples    : qxk
                     "eta",         # etaval               : nxk
                     "lambda",      # loadSamples          : pxk
                     "sigmacol"     # sigmacol (1/sigma^2) : p
  )
  if("all" %in% output) {
    output <- valid_outputs
  } else {
    output <- intersect(output, valid_outputs)
  }
  # -------------------------------------------------------------------------- #
  if(!is.null(W)) {
    if("mu" %in% output) MU <- list()
    if("bmu" %in% output) BMU <- list()
  }
  if("beta" %in% output) BETA <- list()
  if("eta" %in% output) ETA <- list()
  if("lambda" %in% output) LAMBDA <- list()
  if("sigmacol" %in% output) sig <- list()
  K <- rep(NA, sp)
  # -------------------------------------------------------------------------- #
  # ADAPTIVE GIBBS SAMPLING
  # -------------------------------------------------------------------------- #
  ind <- 1
  # start time
  t0 <- proc.time()
  # algorithm starts here
  for (i in 1:nrun) {
    if(verbose == TRUE && i %% 50 == 0) cat(i, ":", k, "active factors\n")
    # ------------------------------------------------------------------------ #
    # 1 - update Z
    if(is.null(W)) {
      Zmean <- tcrossprod(eta, Lambda)  # eta * t(Lambda)
    } else {
      Zmean <- tcrossprod(eta, Lambda) + tcrossprod(W, mu)
    }
    n_unif <- matrix(runif(n * p), nrow = n, ncol = p)
    Z <- truncnorm_lg(y_lower = log(a_y), y_upper = log(a_yp1),
                      mu = Zmean, sigma = sqrt(1 / ps), u_rand = n_unif)
    if(!is.null(W)) {
      # ---------------------------------------------------------------------- #
      # 2 - update b_mu
      b_mu <- update_b_mu_R(X_mean, prec_mu, prec_b, mu, q_mean, c)
      # ---------------------------------------------------------------------- #
      # 3 - update mu
      Z_res <- Z - tcrossprod(eta, Lambda)
      if(c > 1) {
        Qbet <- base::diag(prec_mu, c) + crossprod(W)
      } else {
        Qbet <- prec_mu + crossprod(W)
      }
      for(j in 1:p) {
        mu[j, ] <- update_mu_R(j, Qbet, W, Z_res, ps, b_mu, X_mean, c)
      }
      Z <- Z - tcrossprod(W, mu)
      # ---------------------------------------------------------------------- #
    }
    # 4 - update eta
    eta <- update_eta_R(Lambda, ps, k, Z, n)
    # ------------------------------------------------------------------------ #
    # 5 - update Sigma
    Z_res <- Z - tcrossprod(eta, Lambda)
    ps <- rgamma(p, as + 0.5 * n, bs + 0.5 * colSums(Z_res ^ 2))
    # ------------------------------------------------------------------------ #
    # 6 - update beta
    pred <- X_cov %*% Beta
    logit <- plogis(pred)
    # 6.2 Update phi_L
    Phi_L <- matrix(1, nrow = p, ncol = k)
    logit_phi0 <- logit[which(Phi == 0)]
    which_zero <- which(runif(length(logit_phi0)) < ((1 - logit_phi0) / (1 - logit_phi0 * p_constant)))
    Phi_L[which(Phi == 0)[which_zero]] <- 0
    # 6.3 Polya gamma
    Dt <- matrix(pgdraw::pgdraw(1, pred), nrow = p, ncol = k)
    # 6.4 Update beta_h
    Bh_1 <- diag(rep(sd_beta^2, q_cov) ^ {-1})
    for(h in 1:k) {
      Beta[, h] <- update_beta_R(h, X_cov, Dt, Bh_1, Phi_L, q_cov)
    }
    # ------------------------------------------------------------------------ #
    # 7 - update Lambda and Lambda_star
    etarho <- t(eta) * rho
    for (j in 1:p) {
      Lambda_star[j, ] <- update_Lambda_star_R(j, etarho, Phi, Plam, ps, Z, k)
    }
    Lambda <- t(t(Lambda_star) * sqrt(rho)) * Phi
    # ------------------------------------------------------------------------ #
    # 8.1 - update d
    sdy <- matrix(rep(sqrt(1/ps), n), n, p, byrow = TRUE)
    etalambdastar <- eta[rep(1:n, p), ] * Lambda_star[rep(1:p, each=n), ]
    for(h in 1:k) {
      d[h] <- update_d_R(h, Phi, p, n, rho, etalambdastar, Z, sdy, k, w)
    }
    rho <- rep(1, k)
    rho[d <= seq(1, k)] <- 0
    # 8.2
    Plam <- diag(rgamma(k, a_theta + 0.5 * p, b_theta + 0.5 * colSums(Lambda_star^2)))
    # 8.3
    for(h in 1:(k-1)) {
      v[h] <- rbeta(1, shape1 = 1 + sum(d == h), shape2 = alpha + sum(d > h))
    }
    v[k] <- 1
    w <- v * c(1, cumprod(1 - v[-k]))
    # ------------------------------------------------------------------------ #
    # 9 - update Phi
    pred <- X_cov %*% Beta
    logit <- plogis(pred)
    Phi <- update_Phi_R(rho, logit, p_constant, p, n, eta, Lambda_star, Phi, Z, sdy)
    # ------------------------------------------------------------------------ #
    # save sampled values (after burn-in period)
    if((i %% thin == 0) & (i > burn)) {
      if(!is.null(W)) {
        if("mu" %in% output) MU[[ind]] <- mu
        if("bmu" %in% output) BMU[[ind]] <- b_mu
      }
      if("beta" %in% output) BETA[[ind]] <- Beta
      if("eta" %in% output) ETA[[ind]] <- eta
      if("lambda" %in% output) LAMBDA[[ind]] <- Lambda
      if("sigmacol" %in% output) sig[[ind]] <- ps
      K[ind] <- kstar
      ind <- ind + 1
    }
    # ------------------------------------------------------------------------ #
    # Adaptation
    if((uu[i] < prob[i]) & (i > start_adapt)) {
      active <- which(d > seq(1, k))
      kstar <- length(active)
      if (kstar < k - 1) {
        # set truncation to kstar and subset all variables, keeping only active columns
        k <- kstar + 1
        eta <- cbind(eta[, active, drop = FALSE], rnorm(n))
        vartheta_k <- rgamma(1, a_theta, b_theta)
        Plam <- diag(c(diag(Plam)[active, drop = FALSE], vartheta_k))
        Lambda_star <- cbind(Lambda_star[, active, drop = FALSE], rnorm(p, 0, sd = sqrt(vartheta_k)))
        Phi <- cbind(Phi[, active, drop = FALSE], rbinom(p, size = 1, prob = p_constant))
        rho <- c(rho[active, drop = FALSE], 1)
        Lambda <- cbind(Lambda[, active, drop = FALSE], Lambda_star[, k] * sqrt(rho[k]) * Phi[, k])
        Beta <- cbind(Beta[, active, drop = FALSE], rnorm(q_cov, 0, sd = sqrt(sd_beta)))
        w <- c(w[active, drop = FALSE], 1 - sum(w[active, drop = FALSE]))
        v <- c(v[active, drop = FALSE], 1)  # just to allocate memory
        d <- c(d[active, drop = FALSE], k)
      } else if (k < kmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        k <- k + 1
        eta <- cbind(eta, rnorm(n))
        vartheta_k <- rgamma(1, a_theta, b_theta)
        Plam <- diag(c(diag(Plam), vartheta_k) )
        Lambda_star <- cbind(Lambda_star, rnorm(p, 0, sd = sqrt(vartheta_k)))
        Phi <- cbind(Phi, rbinom(p, size = 1, prob = p_constant))
        rho <- c(rho, 1)
        Lambda <- cbind(Lambda, Lambda_star[, k] * sqrt(rho[k]) * Phi[, k])
        Beta <- cbind(Beta, rnorm(q_cov, 0, sd = sqrt(sd_beta)))
        v[k-1] <- rbeta(1, shape1 = 1, shape2 = alpha)
        v <- c(v, 1)
        w <- v * c(1, cumprod(1 - v[-k]))
        d <- c(d, k)
      }
    }
  }
  # -------------------------------------------------------------------------- #
  out <- lapply(output, function(x) {
    if(!is.null(W)) {
      if(x == "mu") return(MU)          # coefSamples        : pxc
      if(x == "bmu") return(BMU)        # bmuCoefSamples     : qxc
    }
    if(x == "beta") return(BETA)        # loadSamples        : pxk
    if(x == "eta") return(ETA)          # shrinkCoefSamples  : qxk
    if(x == "lambda") return(LAMBDA)    # etaval             : nxk
    if(x == "sigmacol") return(sig)     # sigmacol           : p
  })
  names(out) <- output
  out[["numFactors"]] <- K
  out[["time"]] <- (proc.time() - t0)[1]
  out[["model_prior"]] <- "SIS"
  out[["Y"]] <- Y                       # data               : nxp
  if(!is.null(W)) out[["W"]] <- W       # covariates         : nxc
  out[["X_mean"]] <- X_mean             # metacovariates     : pxq_mean
  out[["X_cov"]]  <- X_cov              # metacovariates     : pxq_cov
  out[["hyperparameters"]] <- list(alpha = alpha, a_theta = a_theta,
                                   b_theta = b_theta, sd_b = sd_b,
                                   sd_mu = sd_mu, sd_beta = sd_beta,
                                   as = as, bs = bs, p_constant = p_constant,
                                   WFormula = WFormula, XmeanFormula = XmeanFormula,
                                   XcovFormula = XcovFormula, y_max = y_max)
  # -------------------------------------------------------------------------- #
  return(out)
  # -------------------------------------------------------------------------- #
}
