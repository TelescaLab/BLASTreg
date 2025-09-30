#' Bayesian Transfer Learning (Oracle, Gaussian) with Horseshoe Shrinkage
#'
#' Oracle BLAST-style regression for Gaussian outcomes with horseshoe prior.
#' Known-informative-set setting; handles target-only edge case.
#' Supports optional empirical Bayes step.
#'
#' @param X Design matrix (target stacked over sources).
#' @param y Outcome vector (target stacked over sources).
#' @param n.vec Integer vector: sizes c(n0, n1, ..., nK). n0 = target size.
#' @param burn Integer burn-in.
#' @param iter Integer saved draws after burn (plus optional EB iterations).
#' @param a,b,s,w Horseshoe / sampler hyperparams (as in your helpers).
#' @param tau Global scale (initial); xi = tau^{-2}.
#' @param sigma2 Error variance (Gaussian).
#' @param alpha Credible interval level.
#' @param iterEBstep Integer EB iterations (0 disables EB).
#' @param sir Logical; passed to EB helper.
#' @param progress_every Print progress every this many iters (0 = silent).
#'
#' @return List with posterior summaries and samples.
#' @export
blast_oracle <- function(
    X, y, n.vec,
    burn = 1000, iter = 5000,
    a = 1/5, b = 10, s = 0.8,
    tau = 1, sigma2 = 1, w = 1,
    alpha = 0.05, iterEBstep = 0, sir = TRUE,
    progress_every = 100
){
  # Basic sizes & splits
  p  <- ncol(X)
  N  <- sum(n.vec)
  n0 <- n.vec[1]
  if (n0 <= 0) stop("n.vec[1] (target size) must be > 0.")
  X0 <- X[seq_len(n0), , drop = FALSE]
  y0 <- y[seq_len(n0)]

  # Source (informative) block
  N_k <- if (length(n.vec) > 1) sum(n.vec[-1]) else 0
  has_sources <- N_k > 0
  if (has_sources) {
    XA <- X[(n0 + 1):N, , drop = FALSE]
    yA <- y[(n0 + 1):N]
  } else {
    XA <- matrix(0, 0, p)
    yA <- numeric(0)
  }

  # Horseshoe / shrinkage initial values
  eta_0    <- rep(1, p)   # local scale^(-2)
  xi_0     <- tau^(-2)    # global scale^(-2)
  sigma2_0 <- sigma2
  kappa_0  <- -0.5
  delta_new <- rep(0, p)

  eta_k    <- rep(1, p)
  xi_k     <- tau^(-2)
  sigma2_k <- sigma2
  kappa_k  <- -0.5

  EB  <- iterEBstep != 0
  nmc <- burn + iterEBstep + iter

  # Precompute crossproducts
  Q0_gauss <- crossprod(X0)
  Qk_gauss <- if (has_sources) crossprod(XA) else NULL
  Q_gauss  <- crossprod(X)

  # Storage
  betaout     <- matrix(0, nrow = nmc, ncol = p)
  wk_out      <- matrix(0, nrow = nmc, ncol = p)
  etaout      <- matrix(0, nrow = nmc, ncol = p)
  eta_kout    <- matrix(0, nrow = nmc, ncol = p)
  xi_0out     <- numeric(nmc)
  xi_kout     <- numeric(nmc)
  sigma2_0out <- numeric(nmc)
  sigma2_kout <- numeric(nmc)
  kappa_0out  <- numeric(nmc)
  kappa_kout  <- numeric(nmc)

  eb_counter <- 0
  sigma2_0_eb_est <- NULL
  sigma2_k_eb_est <- NULL

  ## ---------- TARGET-ONLY SHORT-CIRCUIT ----------
  if (!has_sources) {
    message("BLAST (Oracle, Gaussian): No source studies provided. Proceeding with target-only regression.")

    for (i in seq_len(nmc)) {
      samp <- exact_horseshoe_gibbs_step(
        X = X0, y = y0, eta = eta_0, xi = xi_0, Q = Q0_gauss,
        w = w, s = s, p = p, a = a, b = b,
        sigma2 = sigma2_0, kappa = kappa_0, EB = EB
      )

      beta_new <- samp$new_beta
      eta_0    <- samp$new_eta
      xi_0     <- samp$new_xi
      sigma2_0 <- samp$new_sigma2

      # Save
      betaout[i, ] <- beta_new
      etaout[i, ]  <- eta_0
      xi_0out[i]   <- xi_0
      sigma2_0out[i] <- sigma2_0
      kappa_0out[i] <- kappa_0

      if (progress_every > 0 && i %% progress_every == 0) {
        message(sprintf("Iteration: %d", i))
      }
    }

    burn_all <- burn + iterEBstep
    keep     <- seq.int(burn_all + 1L, nmc)
    beta_s   <- betaout[keep, , drop = FALSE]

    res <- list(
      BetaHat = colMeans(beta_s),
      LeftCI  = apply(beta_s, 2, stats::quantile, probs = alpha/2),
      RightCI = apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2),
      BetaSamples = beta_s
    )
    return(res)
  }

  ## ---------- MAIN LOOP (TARGET + SOURCES) ----------
  for (i in seq_len(nmc)) {
    # --- Auxiliary/shared coefficients w_k update ---
    wk_samp <- exact_horseshoe_gibbs_step(
      X = X, y = c(y0 - drop(X0 %*% delta_new), yA),
      eta = eta_k, xi = xi_k, Q = Q_gauss,
      w = w, s = s, p = p, a = a, b = b,
      sigma2 = sigma2_k, kappa = kappa_k, EB = EB,
      sigma2_0 = sigma2_0, N_k = N_k
    )
    wk_new  <- wk_samp$new_beta
    eta_k   <- wk_samp$new_eta
    xi_k    <- wk_samp$new_xi
    sigma2_k <- wk_samp$new_sigma2

    # --- Target delta update ---
    delta_samp <- exact_horseshoe_gibbs_step(
      X = X0, y = y0 - drop(X0 %*% wk_new),
      eta = eta_0, xi = xi_0, Q = Q0_gauss,
      w = w, s = s, p = p, a = a, b = b,
      sigma2 = sigma2_0, kappa = kappa_0, EB = EB
    )
    delta_new <- delta_samp$new_beta
    eta_0     <- delta_samp$new_eta
    xi_0      <- delta_samp$new_xi
    sigma2_0  <- delta_samp$new_sigma2

    beta <- wk_new + delta_new

    # Save draws
    betaout[i, ]  <- beta
    wk_out[i, ]   <- wk_new
    etaout[i, ]   <- eta_0
    eta_kout[i, ] <- eta_k
    xi_0out[i]    <- xi_0
    xi_kout[i]    <- xi_k
    sigma2_0out[i] <- sigma2_0
    sigma2_kout[i] <- sigma2_k
    kappa_0out[i] <- kappa_0
    kappa_kout[i] <- kappa_k
  }

  # Summaries
  burn_all <- burn + iterEBstep
  keep     <- seq.int(burn_all + 1L, nmc)

  beta_s   <- betaout[keep, , drop = FALSE]
  wk_s     <- wk_out[keep, , drop = FALSE]
  eta_s    <- etaout[keep, , drop = FALSE]
  eta_k_s  <- eta_kout[keep, , drop = FALSE]
  xi0_s    <- xi_0out[keep]
  xik_s    <- xi_kout[keep]
  sig0_s   <- sigma2_0out[keep]
  sigk_s   <- sigma2_kout[keep]
  kappa0_s <- kappa_0out[keep]
  kappak_s <- kappa_kout[keep]

  res <- list(
    BetaHat = colMeans(beta_s),
    LeftCI  = apply(beta_s, 2, stats::quantile, probs = alpha/2),
    RightCI = apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2),
    BetaSamples   = beta_s,
    wkSamples     = wk_s,
    LambdaSamples = 1/sqrt(eta_s),
    TauSamples    = 1/sqrt(xi0_s),
    Sigma2Samples = sig0_s,
    kappa_0 = kappa0_s,
    kappa_k = kappak_s
  )

  return(res)
}
