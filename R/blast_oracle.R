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
#' @param enforce_delta_stronger_shrinkage Logical; if TRUE, during the delta
#'   update the proposal for `log(xi_delta)` is truncated below at `log(xi_wk)`,
#'   ensuring `xi_delta >= xi_wk` and thus stronger/equal shrinkage on contrasts.
#' @param progress_every Integer; print progress every this many iters (0=off).
#'
#'
#'
#'
#'
#' @return List with posterior summaries and samples.
#'
#'
#' #' @examples
#' \dontrun{
#' set.seed(1)
#'
#' # --- Simulate data where all auxiliaries are informative ---
#' M <- 5
#' n.vec <- c(100, rep(100, M))   # target first, then auxiliaries
#' p <- 50; s <- 3
#' size.A0 <- M                   # all auxiliaries informative
#'
#' df <- simulate_multistudy_regression(
#'   p = p, s = s, M = M, size.A0 = size.A0, n.vec = n.vec,
#'   sig.beta = 0.5, sig.delta1 = 0.3, sig.delta2 = 1,
#'   contam_pct = 0.01, type = "gaussian"
#' )
#'
#' # --- Construct oracle gamma (all auxiliaries informative) ---
#' K <- length(n.vec) - 1
#' gamma_oracle <- rep(1L, K)
#'
#' # --- Fit Oracle BLAST ---
#' fit <- blast_oracle(
#'   X = df$X, y = df$y, n.vec = df$n.vec,
#'   gamma = gamma_oracle,
#'   burn = 100, iter = 200,  # keep small for quick examples
#'   a = 1/5, b = 10, s = 0.8,
#'   tau = 1, sigma2 = 1, w = 1, alpha = 0.05
#' )
#'
#' # Inspect estimated coefficients
#' str(fit$BetaHat)
#' }
#'
#' @export
blast_oracle <- function(
  X, y, n.vec,
  burn = 1000, iter = 5000,
  a = 1/5, b = 10, s = 0.8,
  tau = 1, sigma2 = 1, w = 1,
  alpha = 0.05, iterEBstep = 0,
  enforce_delta_stronger_shrinkage = FALSE,
  progress_every = 0
){
  ## --- Dimensions and splits
  p  <- ncol(X)
  N  <- sum(n.vec)
  n0 <- n.vec[1]
  if (n0 <= 0) stop("n.vec[1] (target size) must be > 0.")

  X0 <- X[seq_len(n0), , drop = FALSE]
  y0 <- y[seq_len(n0)]

  N_k <- if (length(n.vec) > 1) sum(n.vec[-1]) else 0
  has_sources <- N_k > 0
  if (has_sources) {
    XA <- X[(n0 + 1):N, , drop = FALSE]
    yA <- y[(n0 + 1):N]
  } else {
    XA <- matrix(0, 0, p)
    yA <- numeric(0)
  }

  ## --- Horseshoe states
  # Target (delta)
  eta_0    <- rep(1, p)
  xi_0     <- tau^(-2)
  sigma2_0 <- sigma2
  kappa_0  <- -0.5
  delta_new <- rep(0, p)

  # Shared/source (wk)
  eta_k    <- rep(1, p)
  xi_k     <- tau^(-2)
  sigma2_k <- sigma2
  kappa_k  <- -0.5

  ## --- Precomputations
  Q0 <- crossprod(X0)
  Q  <- crossprod(X)

  ## --- Storage
  nmc <- burn + iterEBstep + iter
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

  EB <- iterEBstep != 0

  ## --- Target-only short-circuit
  if (!has_sources) {
    message("BLAST (Oracle, Gaussian): No source studies provided. Target-only regression.")
    for (i in seq_len(nmc)) {
      samp <- exact_horseshoe_gibbs_step(
        X = X0, y = y0,
        eta = eta_0, xi = xi_0, Q = Q,
        w = w, s = s, a = a, b = b, p = p,
        sigma2 = sigma2_0, kappa = kappa_0, EB = EB,
        truncated_global_update = FALSE
      )
      beta_new <- samp$new_beta
      eta_0    <- samp$new_eta
      xi_0     <- samp$new_xi
      sigma2_0 <- samp$new_sigma2

      betaout[i, ]   <- beta_new
      etaout[i, ]    <- eta_0
      xi_0out[i]     <- xi_0
      sigma2_0out[i] <- sigma2_0
      kappa_0out[i]  <- kappa_0

      if (progress_every > 0 && i %% progress_every == 0) {
        message(sprintf("Iteration: %d", i))
      }
    }

    keep <- seq.int(burn + iterEBstep + 1L, nmc)
    beta_s <- betaout[keep, , drop = FALSE]
    return(list(
      BetaHat = colMeans(beta_s),
      BetaHatMedian = apply(beta_s, 2, median),
      LeftCI = apply(beta_s, 2, stats::quantile, probs = alpha/2),
      RightCI = apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2),
      Sigma2Hat = mean(sigma2_0out[keep]),
      LambdaHat = colMeans(1/sqrt(etaout[keep, , drop = FALSE])),
      TauHat    = mean(1/sqrt(xi_0out[keep])),
      BetaSamples   = beta_s,
      LambdaSamples = 1/sqrt(etaout[keep, , drop = FALSE]),
      TauSamples    = 1/sqrt(xi_0out[keep]),
      Sigma2Samples = sigma2_0out[keep],
      kappa_0 = kappa_0out[keep]
    ))
  }

  ## --- Main loop (target + sources)
  for (i in seq_len(nmc)) {

    ## (1) Shared/source coefficients wk | delta
    y_new <- c(y0 - drop(X0 %*% delta_new), yA)
    wk_samp <- exact_horseshoe_gibbs_step(
      X = X, y = y_new,
      eta = eta_k, xi = xi_k, Q = Q,
      w = w, s = s, a = a, b = b, p = p,
      sigma2 = sigma2_k, kappa = kappa_k, EB = EB,
      truncated_global_update = FALSE
    )
    wk_new   <- wk_samp$new_beta
    eta_k    <- wk_samp$new_eta
    xi_k     <- wk_samp$new_xi
    sigma2_k <- wk_samp$new_sigma2

    ## (2) Target delta | wk_new
    # If enforcing stronger shrinkage on contrasts: truncate log-proposal at log(xi_k)
    delta_samp <- exact_horseshoe_gibbs_step(
      X = X0, y = y0 - drop(X0 %*% wk_new),
      eta = eta_0, xi = xi_0, Q = Q0,
      w = w, s = s, a = a, b = b, p = p,
      sigma2 = sigma2_0, kappa = kappa_0, EB = EB,
      truncated_global_update = enforce_delta_stronger_shrinkage,
      truncation_threshold    = if (enforce_delta_stronger_shrinkage) xi_k else NULL
    )
    delta_new <- delta_samp$new_beta
    eta_0     <- delta_samp$new_eta
    xi_0      <- delta_samp$new_xi
    sigma2_0  <- delta_samp$new_sigma2

    beta <- wk_new + delta_new

    ## Save draws
    betaout[i, ]    <- beta
    wk_out[i, ]     <- wk_new
    etaout[i, ]     <- eta_0
    eta_kout[i, ]   <- eta_k
    xi_0out[i]      <- xi_0
    xi_kout[i]      <- xi_k
    sigma2_0out[i]  <- sigma2_0
    sigma2_kout[i]  <- sigma2_k
    kappa_0out[i]   <- kappa_0
    kappa_kout[i]   <- kappa_k

    ## (3) EB update *inside* the sampler (oracle: no marginal adjustment)
    if (EB && i == burn + iterEBstep) {
      tmp_k <- eb_em_max(
        kappa = kappa_k, N = N_k,
        xi_out = xi_kout, sigma2_out = sigma2_kout,
        iterEBstep = iterEBstep, eb_counter = 0, i = i
      )
      kappa_k <- tmp_k$kappa

      tmp_0 <- eb_em_max(
        kappa = kappa_0, N = n0,
        xi_out = xi_0out, sigma2_out = sigma2_0out,
        iterEBstep = iterEBstep, eb_counter = 0, i = i
      )
      kappa_0 <- tmp_0$kappa

      if (progress_every > 0) {
        message(sprintf("EB updated: kappa_0=%.3f, kappa_k=%.3f", kappa_0, kappa_k))
      }
    }

    if (progress_every > 0 && i %% progress_every == 0) {
      message(sprintf("Iteration: %d", i))
    }
  }

  ## --- Summaries
  keep <- seq.int(burn + iterEBstep + 1L, nmc)

  beta_s   <- betaout[keep, , drop = FALSE]
  wk_s     <- wk_out[keep, , drop = FALSE]
  eta_s    <- etaout[keep, , drop = FALSE]
  xi0_s    <- xi_0out[keep]
  sig0_s   <- sigma2_0out[keep]

  res <- list(
    BetaHat = colMeans(beta_s),
    BetaHatMedian = apply(beta_s, 2, median),
    LeftCI  = apply(beta_s, 2, stats::quantile, probs = alpha/2),
    RightCI = apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2),
    Sigma2Hat = mean(sig0_s),
    LambdaHat = colMeans(1/sqrt(eta_s)),
    TauHat    = mean(1/sqrt(xi0_s)),
    BetaSamples   = beta_s,
    wkSamples     = wk_s,
    LambdaSamples = 1/sqrt(eta_s),
    TauSamples    = 1/sqrt(xi0_s),
    Sigma2Samples = sig0_s,
    kappa_0 = kappa_0out[keep],
    kappa_k = kappa_kout[keep]
  )

  return(res)
}
