# ---------------------------------------------------------------------------- #

#' Update b_mu in the Adaptive Gibbs Sampler
#'
#' @param X A pxq matrix.
#' @param prec_mu A positive number.
#' @param prec_b A positive number.
#' @param mu A pxc matrix.
#' @param q An integer.
#' @param c An integer.
#'
#' @return A qxc matrix.
#'
#' @importFrom stats rnorm
update_b_mu_R <- function(X, prec_mu, prec_b, mu, q, c) {
  Trsg <- X * prec_mu
  if(q > 1) {
    Vgmbeta1 <- base::diag(prec_b, q) + crossprod(Trsg, X)
  } else {
    Vgmbeta1 <- prec_b + crossprod(Trsg, X)
  }
  eigs <- eigen(Vgmbeta1)
  if(all(eigs$values > 1e-6)) {
    Tmat <- sqrt(eigs$values) * t(eigs$vectors)
  } else {
    Tmat <- chol(Vgmbeta1)
  }
  R <- qr.R(qr(Tmat))
  S <- solve(R)
  Vgmbeta <- S %*% t(S)
  Mgmbeta <- crossprod(mu, Trsg %*% Vgmbeta)  # cxq
  b_mu <- t(Mgmbeta + matrix(rnorm(c * q), nrow = c, ncol = q) %*% t(S))  # qxc
  b_mu
}

# ---------------------------------------------------------------------------- #

#' Update the hth row of mu in the Adaptive Gibbs Sampler
#'
#' @param j DA SCRIVERE
#' @param Qbet DA SCRIVERE
#' @param W DA SCRIVERE
#' @param Z_res DA SCRIVERE
#' @param ps DA SCRIVERE
#' @param b_mu DA SCRIVERE
#' @param Xcov DA SCRIVERE
#' @param c DA SCRIVERE
#'
#' @return A 1xc matrix.
#'
#' @importFrom stats rnorm
update_mu_R <- function(j, Qbet, W, Z_res, ps, b_mu, Xcov, c) {
  Lbet <- t(chol(Qbet))
  bbet <- t(W) %*% Z_res[, j] + ps[j] * crossprod(b_mu, t(Xcov)[, j])
  # mean
  vbet <- forwardsolve(Lbet, bbet)
  mbet <- backsolve(t(Lbet), vbet)
  # var
  zbet <- rnorm(c)
  ybet <- backsolve(t(Lbet), zbet)
  t(ybet + mbet)
}

# ---------------------------------------------------------------------------- #

#' Update eta in the Adaptive Gibbs Sampler
#'
#' @param Lambda A pxk matrix.
#' @param ps DA SCRIVERE
#' @param k Number of columns of eta.
#' @param Z A nxp matrix.
#' @param n Number of rows of eta.
#'
#' @return A nxk matrix.
#'
#' @importFrom stats rnorm
update_eta_R <- function(Lambda, ps, k, Z, n) {
  Lmsg <- Lambda * ps
  Veta1 <- base::diag(k) + crossprod(Lmsg, Lambda)
  eigs <- eigen(Veta1)
  if(all(eigs$values > 1e-6)) {
    Tmat <- sqrt(eigs$values) * t(eigs$vectors)
  } else {
    Tmat <- chol(Veta1)
  }
  R <- qr.R(qr(Tmat))
  S <- solve(R)
  Veta <- S %*% t(S)
  Meta <- Z %*% Lmsg %*% Veta
  eta <- Meta + matrix(rnorm(n * k), nrow = n, ncol = k) %*% t(S)
  eta
}

# ---------------------------------------------------------------------------- #

#' Update the hth column of beta in the Adaptive Gibbs Sampler
#'
#' @param h DA SCRIVERE
#' @param Xcov DA SCRIVERE
#' @param Dt DA SCRIVERE
#' @param Bh_1 DA SCRIVERE
#' @param Phi_L DA SCRIVERE
#' @param q DA SCRIVERE
#'
#' @return A qx1 matrix.
#'
#' @importFrom stats rnorm
update_beta_R <- function(h, Xcov, Dt, Bh_1, Phi_L, q) {
  Qbeta <- t(Dt[, h] * Xcov) %*% Xcov + Bh_1
  bbeta <- crossprod(Xcov, (Phi_L[, h] - 0.5))
  Lbeta <- t(chol(Qbeta))
  # mean
  vbeta <- forwardsolve(Lbeta, bbeta)
  mbeta <- backsolve(t(Lbeta), vbeta)
  # var
  zbeta <- rnorm(q)
  ybeta <- backsolve(t(Lbeta), zbeta)
  t(ybeta + mbeta)
}

# ---------------------------------------------------------------------------- #

#' Update the jth row of Lambda_star in the Adaptive Gibbs Sampler
#'
#' @param j DA SCRIVERE
#' @param etarho DA SCRIVERE
#' @param Phi DA SCRIVERE
#' @param Plam DA SCRIVERE
#' @param ps DA SCRIVERE
#' @param Z DA SCRIVERE
#' @param k DA SCRIVERE
#'
#' @return A 1xk matrix.
#'
#' @importFrom stats rnorm
update_Lambda_star_R <- function(j, etarho, Phi, Plam, ps, Z, k) {
  etaj <- t(etarho * Phi[j,])
  eta2 <- crossprod(etaj)
  Qlam <- Plam + ps[j] * eta2
  blam <- ps[j] * crossprod(etaj, Z[,j])
  Llam <- t(chol(Qlam))
  # mean
  vlam <- forwardsolve(Llam, blam)
  mlam <- backsolve(t(Llam), vlam)
  # var
  zlam <- rnorm(k)
  ylam <- backsolve(t(Llam), zlam)
  t(ylam + mlam)
}

# ---------------------------------------------------------------------------- #

#' Update the hth element of d in the Adaptive Gibbs Sampler
#'
#' @param h DA SCRIVERE
#' @param Phi DA SCRIVERE
#' @param p DA SCRIVERE
#' @param n DA SCRIVERE
#' @param rho DA SCRIVERE
#' @param etalambdastar DA SCRIVERE
#' @param Z DA SCRIVERE
#' @param sdy DA SCRIVERE
#' @param k DA SCRIVERE
#' @param w DA SCRIVERE
#'
#' @return An integer in 1, ..., k.
#'
#' @importFrom stats dnorm rmultinom
update_d_R <- function(h, Phi, p, n, rho, etalambdastar, Z, sdy, k, w) {
  phirho <- Phi[rep(1:p, each = n), -h] * sqrt(rho[-h])
  muijh0 <- rowSums(etalambdastar[, -h] * phirho)
  lnorm0 <- sum(dnorm(Z, mean = muijh0, sd = sdy, log=T))
  lnorm1 <- sum(dnorm(Z, mean = muijh0 + etalambdastar[,h] * Phi[rep(1:p, each=n), h], sd = sdy, log=T))
  # adjust the scale
  mlnorm <- max(c(lnorm0, lnorm1))
  lnorm0 <- lnorm0 - mlnorm
  lnorm1 <- lnorm1 - mlnorm
  prob_h <- exp(c(rep(lnorm0, h), rep(lnorm1, k - h)) + log(w))
  if(sum(prob_h) == 0) prob_h <- c(rep(0, k - 1), 1)
  which(rmultinom(n=1, size=1, prob=prob_h) == 1)
}

# ---------------------------------------------------------------------------- #

#' Update Phi in the Adaptive Gibbs Sampler
#'
#' @param rho DA SCRIVERE
#' @param logit DA SCRIVERE
#' @param p_constant DA SCRIVERE
#' @param p DA SCRIVERE
#' @param n DA SCRIVERE
#' @param eta DA SCRIVERE
#' @param Lambda_star DA SCRIVERE
#' @param Phi DA SCRIVERE
#' @param Z DA SCRIVERE
#' @param sdy DA SCRIVERE
#'
#' @return A pxk matrix.
#'
#' @importFrom stats dnorm runif
update_Phi_R <- function(rho, logit, p_constant, p, n,
                       eta, Lambda_star, Phi, Z, sdy) {
  # select active factors
  wr  <- which(rho == 1)
  wr0 <- which(rho == 0)
  # update for inactive factors
  for(h in wr0) {
    lp_phi0 <- log(1 - logit[, h] * p_constant)
    lp_phi1 <- log(logit[, h] * p_constant)
    sumlog <- apply(cbind(lp_phi0, lp_phi1), 1, matrixStats::logSumExp)
    Phi[, h] <- round(runif(p) < exp(lp_phi1 - sumlog))
  }
  # update for active factors
  etalambdastar <- eta[rep(1:n, p), wr] * Lambda_star[rep(1:p, each=n), wr]
  ih <- 0
  for(h in wr) {
    ih <- ih + 1
    muijh <- rowSums(etalambdastar[, -ih] * Phi[rep(1:p, each=n), wr[-ih]])
    lnorm0 <- tapply(dnorm(Z, mean = muijh, sd = sdy, log=T), rep(1:p, each=n), FUN = sum)
    lnorm1 <- tapply(dnorm(Z, mean = muijh + etalambdastar[, ih] * rho[h], sd = sdy, log=T), rep(1:p, each=n), FUN = sum)
    # adjust scale
    mlnorm <- max(c(lnorm0, lnorm1))
    lnorm0 <- lnorm0 - mlnorm
    lnorm1 <- lnorm1 - mlnorm
    lp_phi0 <- lnorm0 + log(1 - logit[, h] * p_constant)
    lp_phi1 <- lnorm1 + log(logit[, h] * p_constant)
    sumlog <- apply(cbind(lp_phi0, lp_phi1), 1, matrixStats::logSumExp)
    Phi[, h] <- round(runif(p) < exp(lp_phi1 - sumlog))
  }
  Phi
}

# ---------------------------------------------------------------------------- #
