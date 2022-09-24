#' AGS for GIF models with SIS prior (C++)
#'
#' @description
#' Implementation in C++ of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
#'
#' @details
#' Suppose an \eqn{n\times p} matrix \eqn{y} of counts is available.
#' We consider a count-valued stochastic process \eqn{y_{ij}:\mathbb{W}\to\mathbb{N}}, where \eqn{\mathbb{W}} is the sample space and \eqn{\mathbb{N}=\{0,\ldots,\infty\}}, \eqn{i=1,\ldots,n} and \eqn{j=1,\ldots,p}.
#' We introduce continuous-valued process \eqn{y^*_{ij}:\mathbb{X}\to\mathbb{T}}, \eqn{\mathbb{T}\subseteq\mathbb{R}} related to the observed count-valued data \eqn{y_{ij}} via
#' \deqn{y_{ij} = h(y^*_{ij})}
#' where \eqn{h:\mathbb{T}\to\mathbb{N}} is a rounding operator that sets \eqn{y_{ij}(w)=t} when \eqn{y^*_{ij}(w)\in\mathbb{A}_t} and \eqn{\{\mathbb{A}_t\}^\infty_{t=1}} is a known partition of \eqn{\mathbb{T}}.
#' We introduce latent variables \eqn{z_{ij}}, defined as \eqn{\log(y^*_{ij}) = z_{ij}} and modeled through an Infinite Factor model as follows
#' \deqn{z_{ij} = w_i^\top\mu_j+\epsilon_{ij},}
#' where \eqn{w_i\in\mathbb{R}^c} are the covariates of the \eqn{i}th observation, \eqn{\mu_j\in\mathbb{R}^c} quantifies the effect of the covariates on the \eqn{j}th column of \eqn{y}, and
#' \deqn{\epsilon_i=(\epsilon_{i1},\ldots,\epsilon_{ip})^\top\sim N_p(0, \Omega)}
#' The matrix \eqn{\Omega = var(\epsilon_i)} can be expressed as
#' \deqn{\Omega = \Lambda\Lambda^\top + \Sigma}
#' where the \eqn{h}th column of \eqn{\Lambda\in\mathbb{R}^{p\times k}} quantifies the effect of the \eqn{h}th latent factor on the columns of \eqn{y} and \eqn{\Sigma=diag(\sigma^2_1,\ldots,\sigma^2_p)} with \eqn{\sigma^{-2}_j\sim Ga(a_{\sigma},b_{\sigma})}, \eqn{j=1,\ldots,p}.
#'
#' The focus is on a new class of generalized infinite factor models induced through a novel class of priors for \eqn{\Lambda} that allows infinitely many factors, \eqn{k = \infty}. In particular, we let
#' \deqn{\lambda_{jh}|\theta_{jh}\sim N(0,\theta_{jh}), \quad \theta_{jh}=\tau_0\gamma_h\phi_{jh}, \quad \tau_0\sim f_{\tau_0}, \quad \gamma_h\sim f_{\gamma_h}, \quad \phi_{jh}\sim f_{\phi_j},}
#' where \eqn{f_{\tau_0}}, \eqn{f_{\gamma_h}} and \eqn{f_{\phi_j}} are supported in \eqn{[0,\infty)} with positive probability mass on \eqn{(0,\infty)}. The local \eqn{\phi_{jh}}, column-specific \eqn{\gamma_h}, and global \eqn{\tau_0} scales are all independent a priori.
#'
#' We define a non-exchangeable structure, called \emph{Structured Increasing Shrinkage} (\emph{SIS}) Process, that includes meta-covariates \eqn{x\in\mathbb{R}^{p\times q}} informing the sparsity structure of \eqn{\Lambda}; we specify
#' \deqn{\tau_0=1, \quad\gamma_h=\theta_h\rho_h, \quad\phi_{jh}|\beta_h\sim Ber(logit^{-1}(x_j^\top\beta_h)c_p),}
#' \deqn{\theta_h^{-1}\sim Ga(a_\theta,b_\theta), \quad a_\theta>1, \quad\rho_h=Ber(1-\pi_h), \quad \beta_h\sim N_q(0, \sigma^2_\beta I_q),}
#' where we assume the link \eqn{g(x)=logit^{-1}(x)c_p}, with \eqn{logit^{-1}(x)=e^x/(1+e^x)} and \eqn{c_p\in(0,1)} a possible offset.
#' The parameter \eqn{\pi_h=pr(\lambda_h=0)} follows a stick-breaking construction,
#' \deqn{\pi_h=\sum_{l=1}^hw_l, \quad w_l=v_l\prod_{m=1}^{l-1}(1-v_m), \quad v_m\sim Be(1,\alpha),}
#'
#' Posterior inference is conducted via Markov chain Monte Carlo sampling.
#' Following common practice in infinite factor models, we use an Adaptive Gibbs Sampler, which attempts to infer the best truncation level \eqn{H} while drawing from the posterior distribution of the parameters.
#' The value of \eqn{H} is adapted only at some Gibbs iterations by discarding redundant factors and, if no redundant factors are identified, by adding a new factor by sampling its parameters from the prior distribution.
#' The probability of occurrence of an adaptive iteration \eqn{t} as equal to \eqn{p(t)=\exp(b_0+b_1t)}, where \eqn{b_0} and \eqn{b_1} are positive constants, such that frequency of adaptation decreases.
#'
#' @param Y A \eqn{n\times p} matrix \eqn{y} of counts.
#' @param X A matrix \eqn{x} of meta-covariates having \eqn{p} rows; the variables must be numeric or factors.
#' @param W A matrix \eqn{w} of covariates having \eqn{n} rows; the variables must be numeric or factors.
#' @param seed Seed. Default is 28.
#' @param stdx Logical: if \code{TRUE}, numeric meta-covariates are standardized; by default meta-covariates are standardized.
#' @param stdw Logical: if \code{TRUE}, numeric covariates are standardized; by default covariates are standardized.
#' @param XFormula Formula specifying  the meta-covariates inserted in the model; by default all are considered.
#' @param WFormula Formula specifying  the covariates inserted in the model; by default all are considered.
#' @param kinit An integer minimun number of latent factors. Default is \code{min(floor(log(p)*kval), p)}.
#' @param kmax Maximum number of latent factors. Default is \code{p+1}.
#' @param kval An integer number used to calculate the default value of \code{kinit}. Default is 6.
#' @param nrun An integer number of iterations. Default is 100.
#' @param burn An integer number of burning iterations (number of iterations to discard). Default is \code{round(nrun/4)}.
#' @param thin An integer thinning value (number of iterations to skip between saving iterations). Default is 1.
#' @param start_adapt An integer number of iterations before adaptation. Default is 50.
#' @param y_max A fixed and known upper bound for the values in \code{Y}. Default is \code{Inf}.
# Parameters for SIS:
#' @param b0,b1 Positive constants for the adaptive probability \eqn{p(t)=\exp(b_0+b_1t)}. Default is \eqn{b_0=1} and \eqn{b_1=5\times 10^{-4}}.
#' @param sd_b Standard deviation for \eqn{b_m}. Default is \eqn{\sigma_b=1}.
#' @param sd_mu Standard deviation for \eqn{\mu_j}. Default is \eqn{\sigma_\mu=1}.
#' @param sd_beta Standard deviation for \eqn{\beta_h}. Default is \eqn{\sigma_\beta=1}.
#' @param a_theta,b_theta Shape (\code{a_theta}) and rate (\code{b_theta}) parameters of the gamma prior distribution of \eqn{\theta_{jh}^{-1}}. Default is \eqn{a_{\theta}=b_{\theta}=1}.
#' @param as,bs Shape (\code{as}) and rate (\code{bs}) parameters of the gamma prior distribution of \eqn{\sigma_j^{-2}}. Default is \eqn{a_{\sigma}=b_{\sigma}=1}.
#' @param p_constant Factor probability constant. Default is \eqn{c_p=10e\log(p)/p}.
#' @param alpha Non-negative parameter of the Beta prior distribution of \eqn{v_m}, \eqn{Be(1,\alpha)}. Default is \eqn{\alpha=5}.
# Other settings:
#' @param output A vector containing the names of the parameters for which you want to save the draws from the posterior distribution. The possible valid strings are \code{"mu"}, \code{"bmu"}, \code{"beta"}, \code{"eta"}, \code{"lambda"} and \code{"sigmacol"}. Default is \code{"all"}, which is equivalent to writing \code{c("mu", "bmu", "beta", "eta", "lambda", "sigmacol")}.
#' @param verbose Logical: if \code{TRUE}, print the number of active factors every 50 iterations. Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{numFactors}: a vector containing the number of active factors at each saved iteration.
#' \item \code{mu}: a list containing draws from the posterior distribution of \eqn{\mu=(\mu_1,\ldots,\mu_p)^\top\in\mathbb{R}^{p\times c}}.
#' \item \code{bmu}: a list containing draws from the posterior distribution of \eqn{b_{\mu}=(b_{\mu,1},\ldots,b_{\mu,q})^\top\in\mathbb{R}^{c\times q}}.
#' \item \code{beta}: a list containing draws from the posterior distribution of \eqn{\beta=(\beta_1,\ldots,\beta_k)\in\mathbb{R}^{q\times k}}.
#' \item \code{eta}: a list containing draws from the posterior distribution of \eqn{\eta=(\eta_1,\ldots,\eta_n)^\top\in\mathbb{R}^{n\times k}}.
#' \item \code{lambda}: a list containing draws from the posterior distribution of \eqn{\Lambda\in\mathbb{R}^{p\times k}}.
#' \item \code{sigmacol}: a list containing draws from the posterior distribution of \eqn{(\sigma_1^{-2},\ldots,\sigma_p^{-2})\in\mathbb{R}^{p}}.
#' \item \code{time}:
#' \item \code{model_prior}: the name of the class of priors for \eqn{\Lambda}.
#' \item \code{Y}: the \eqn{n\times p} matrix of counts provided as input to the function.
#' \item \code{W}: the \eqn{p\times q} matrix of meta-covariates obtained as a result of variable selection, using \code{WFormula}, and conversion of factors into dichotomous variables.
#' \item \code{X}: the \eqn{n\times c} matrix of covariates obtained as a result of variable selection, using \code{XFormula}, and conversion of factors into dichotomous variables.
#' \item \code{hyperparameters}: a list containing the hyperparameters provided as input to the function.
#' }
#'
#' @seealso The function \code{\link{lposterior}} compute the log-posterior probabilities of part or all the MCMC iterations.. Two alternative implementations of the Adaptive Gibbs Sampler are \code{\link{AGS_SIS_R}} and \code{\link{AGS_SIS_RC}}.
#'
#' @importFrom stats formula model.matrix plogis rbeta rbinom rgamma rnorm runif
#'
#' @import RcppArmadillo
#'
#' @export
AGS_SIS <- function(Y,
                    X = NULL, W = NULL,
                    seed = 28,
                    stdx = TRUE, stdw = TRUE,
                    WFormula = formula("~ ."),
                    XFormula = formula("~ ."),
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
  if(!is.data.frame(Y) && !is.matrix(Y)) {
    stop("'Y' not valid: it must be a matrix or a dataframe.")
  }
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  # -------------------------------------------------------------------------- #
  if(is.null(X)) {
    X <- matrix(1, nrow = p, ncol = 1)
  } else {
    if(!is.data.frame(X) && !is.matrix(X)) {
      stop("'X' not valid: it must be a matrix or a dataframe.")
    }
    if(p != nrow(X)) {
      stop("'Y' and 'X' not compatible: the number of columns of 'Y' must be equal to the number of rows of 'X'.")
    }
    if((length(stdx) != 1) || !is.logical(stdx)) {
      stop("'stdx' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdx) {
      is.fact.x <- sapply(X, is.factor)
      X[, is.fact.x == FALSE] <- scale(X[, is.fact.x == FALSE])
      if(is.data.frame(X) & dim(X)[2] > 1) {
        X <- model.matrix(XFormula, X)
      }
    }
  }
  q <- dim(X)[2]
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
    if(stdw) {
      is.fact.w <- sapply(W, is.factor)
      W[, is.fact.w == FALSE] <-  scale(W[, is.fact.w == FALSE])
      if(is.data.frame(W) & dim(W)[2] > 1) {
        W <- model.matrix(WFormula, W)
      }
    }
    c <- ncol(W)
  } else {
    W <- matrix(1, nrow = 1, ncol = 1)
    c <- 1
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
    mu <- matrix(rnorm(c * p, 0, sd_mu), nrow = p, ncol = c)  # mean coeff of the data
    b_mu <- matrix(rnorm(q * c), nrow = c, ncol = q)          # x effects on mu coeff
    # precision of mu and b_mu
    prec_b  <- 1 / (sd_b)  ^ 2
    prec_mu <- 1 / (sd_mu) ^ 2
  } else {
    # not used but necessary to call the C++ function
    mu <- b_mu <- matrix(1, nrow = 1, ncol = 1)
    prec_b <- prec_mu <- 0.1
  }
  # Initialize lambda star (pxq)
  Lambda_star <- matrix(rnorm(p * k), nrow = p, ncol = k)  # loading matrix
  # Initialize eta (nxk)
  eta <- matrix(rnorm(n * k), nrow = n, ncol = k)          # latent factors
  # Initialize Beta (qxk) and pred (pxk)
  Beta <- matrix(rnorm(q * k), nrow = q, ncol = k)  # traits effect on local shrinkage
  pred <- X %*% Beta                                # local shrinkage coefficients
  logit <- plogis(pred)
  # Initialize Phi pxk
  Phi <- matrix(rbinom(p * k, size = 1, prob = p_constant), nrow = p, ncol = k)
  # Initialize pi_h, h = 1, ..., k
  v <- c(rbeta(k - 1, shape1 = 1, shape2 = alpha), 1)
  w <- v * c(1, cumprod(1 - v[-k]))  # product up to  l - 1
  d <- rep(k - 1, k)                 # augmented data
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
  out <- list("numFactors" = NA)
  if(!Wnull) {
    if("mu" %in% output) out["mu"] <- NA
    if("bmu" %in% output) out["bmu"] <- NA
  }
  if("beta" %in% output) out["beta"] <- NA
  if("eta" %in% output) out["eta"] <- NA
  if("lambda" %in% output) out["lambda"] <- NA
  if("sigmacol" %in% output) out["sigmacol"] <- NA
  # -------------------------------------------------------------------------- #
  # start time
  t0 <- proc.time()
  # -------------------------------------------------------------------------- #
  # ADAPTIVE GIBBS SAMPLING
  # -------------------------------------------------------------------------- #
  out <- Rcpp_AGS_SIS(alpha, as, a_y, a_yp1, a_theta,
                      Beta, bs, b0, b1, burn, b_mu, b_theta,
                      c,
                      d,
                      eta,
                      k,  kmax, kstar,
                      Lambda, Lambda_star, logit,
                      mu,
                      n, nrun,
                      out,
                      p, Phi, Plam, prec_b, prec_mu, pred, prob, ps, p_constant,
                      q,
                      rho,
                      sd_b, sd_beta, sd_mu, sp, start_adapt,
                      thin,
                      uu,
                      v, verbose,
                      w, W, Wnull,
                      X)
  # -------------------------------------------------------------------------- #
  for (it in 1:sp) {
    if("sigmacol" %in% output) out$sigmacol[[it]] <- c(out$sigmacol[[it]])
  }
  out["numFactors"] <- c(out["numFactors"])
  out["time"] <- (proc.time() - t0)[1]
  out[["model_prior"]] <- "SIS"
  out[["Y"]] <- Y                       # data               : nxp
  if(!Wnull) out[["W"]] <- W            # covariates         : nxc
  out[["X"]] <- X                       # metacovariates     : pxq
  out[["hyperparameters"]] <- list(alpha = alpha, a_theta = a_theta,
                                   b_theta = b_theta, sd_b = sd_b,
                                   sd_mu = sd_mu, sd_beta = sd_beta,
                                   as = as, bs = bs, p_constant = p_constant,
                                   WFormula = WFormula, XFormula = XFormula,
                                   y_max = y_max)
  # -------------------------------------------------------------------------- #
  return(out)
  # -------------------------------------------------------------------------- #
}
