#' Bayesian Transfer Learning (Oracle) with Horseshoe Shrinkage
#'
#' Unified Gaussian / Binomial (logistic with Pólya–Gamma) oracle BLAST-style regression.
#' Known-informative-set setting; handles target-only edge case. Supports optional EB step.
#'
#' @param X Design matrix (target stacked over sources).
#' @param y Outcome vector (target stacked over sources).
#' @param n.vec Integer vector: sizes c(n0, n1, ..., nK). n0 = target size.
#' @param burn Integer burn-in.
#' @param iter Integer saved draws after burn (plus optional EB iterations).
#' @param a,b,s,w Horseshoe / sampler hyperparams (as in your helpers).
#' @param tau Global scale (initial); xi = tau^{-2}.
#' @param sigma2 Error var (Gaussian only; ignored in logistic).
#' @param alpha Credible interval level.
#' @param iterEBstep Integer EB iterations (0 disables EB).
#' @param sir Logical; passed to EB helper.
#' @param model "gaussian" or "binomial".
#' @param progress_every Print progress every this many iters (0 = silent).
#'
#' @return List with posterior summaries and samples (sigma2 fields NULL in logistic).
#' @export

blast_oracle <- function(
    X, y, n.vec,
    burn = 1000, iter = 5000,
    a = 1/5, b = 10, s = 0.8,
    tau = 1, sigma2 = 1, w = 1,
    alpha = 0.05, iterEBstep = 0, sir = TRUE,
    model = c("gaussian","binomial"),
    progress_every = 100
){
  model <- match.arg(tolower(model), c("gaussian","binomial"))

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

  # Horseshoe / shrinkage initial values (shared across models)
  # Target (delta) block
  eta_0   <- rep(1, p)                 # local scale^(-2), consistent w/ your code
  xi_0    <- tau^(-2)                  # global scale^(-2)
  sigma2_0 <- sigma2                   # Gaussian only; ignored in logistic
  kappa_0  <- -0.5                     # your default for EB (used in both branches)
  delta_new <- rep(0, p)

  # Auxiliary (wk) block
  eta_k   <- rep(1, p)
  xi_k    <- tau^(-2)
  sigma2_k <- sigma2
  kappa_k  <- -0.5

  # Empirical Bayes toggle
  EB <- iterEBstep != 0
  nmc <- burn + iterEBstep + iter

  # Precompute Q's for Gaussian case (constant across iterations)
  Q0_gauss <- if (model == "gaussian") crossprod(X0) else NULL
  Qk_gauss <- if (model == "gaussian" && has_sources) crossprod(XA) else NULL
  Q_gauss  <- if (model == "gaussian") crossprod(X) else NULL

  # Initialize PG omegas for logistic
  if (model == "binomial") {
    # Start with ones (common choice)
    omega_0 <- rep(1, n0)
    omega_A <- if (has_sources) rep(1, N_k) else numeric(0)
  }

  # Storage
  betaout     <- matrix(0, nrow = nmc, ncol = p)
  wk_out      <- matrix(0, nrow = nmc, ncol = p)
  etaout      <- matrix(0, nrow = nmc, ncol = p)       # target local scales
  eta_kout    <- matrix(0, nrow = nmc, ncol = p)       # aux local scales
  xi_0out     <- numeric(nmc)
  xi_kout     <- numeric(nmc)
  sigma2_0out <- if (model == "gaussian") numeric(nmc) else NULL
  sigma2_kout <- if (model == "gaussian") numeric(nmc) else NULL
  kappa_0out  <- numeric(nmc)
  kappa_kout  <- numeric(nmc)

  eb_counter <- 0
  sigma2_0_eb_est <- NULL
  sigma2_k_eb_est <- NULL

  ## ---------- TARGET-ONLY SHORT-CIRCUIT ----------
  if (!has_sources) {
    message("BLAST (Oracle): No source studies provided. Proceeding with target-only regression.")

    for (i in seq_len(nmc)) {
      if (model == "gaussian") {
        # Standard Gaussian horseshoe step on target
        samp <- exact_horseshoe_gibbs_step(
          X = X0, y = y0, eta = eta_0, xi = xi_0, Q = Q0_gauss,
          w = w, s = s, p = p, a = a, b = b,
          sigma2 = sigma2_0, kappa = kappa_0,
          EB = EB
        )
      } else {
        # Logistic target-only: weight X and transform y using PG
        X0_w <- sweep(X0, 1, sqrt(omega_0), `*`)
        y0_w <- (y0 - 1/2) / sqrt(omega_0)  # no wk/delta split here: beta = delta
        Q0   <- crossprod(X0_w)

        samp <- exact_horseshoe_gibbs_step(
          X = X0_w, y = y0_w, eta = eta_0, xi = xi_0, Q = Q0,
          w = w, s = s, p = p, a = a, b = b,
          sigma2 = NULL, kappa = kappa_0,
          EB = EB
        )
      }

      # Extract
      beta_new <- samp$new_beta
      eta_0    <- samp$new_eta
      xi_0     <- samp$new_xi
      if (model == "gaussian") sigma2_0 <- samp$new_sigma2

      # Save
      betaout[i, ] <- beta_new
      etaout[i, ]  <- eta_0
      xi_0out[i]   <- xi_0
      if (model == "gaussian") sigma2_0out[i] <- sigma2_0
      kappa_0out[i] <- kappa_0

      # EB step (after burn + iterEBstep)
      if (EB && i == burn + iterEBstep) {
        tmp_0 <- eb_em_max(
          kappa = kappa_0, N = n0,
          xi_out = xi_0out, sigma2_out = if (model == "gaussian") sigma2_0out else rep(NA_real_, nmc),
          iterEBstep = iterEBstep, eb_counter = eb_counter, i = i, sir = sir
        )
        kappa_0 <- tmp_0$kappa
        sigma2_0_eb_est <- if (model == "gaussian") tmp_0$sigma2_est else NULL
        eb_counter <- eb_counter + 1
        if (progress_every > 0) {
          val1 <- (exp(kappa_0) / (1 + exp(kappa_0))) * p
          message(sprintf("EB updated (target): kappa_0=%.3f, E[sparsity]=%.2f", kappa_0, val1))
        }
      }

      # Update omegas for logistic
      if (model == "binomial") {
        # beta_new IS the target coefficient vector
        # PG parameter is linear predictor x^T beta
        for (j in seq_len(n0)) {
          omega_0[j] <- pgdraw::pgdraw(1, drop(X0[j, , drop = FALSE] %*% beta_new))
        }
      }

      if (progress_every > 0 && i %% progress_every == 0) {
        message(sprintf("Iteration: %d", i))
      }
    }

    # Summaries
    burn_all <- burn + iterEBstep
    keep     <- seq.int(burn_all + 1L, nmc)
    beta_s   <- betaout[keep, , drop = FALSE]
    eta_s    <- etaout[keep, , drop = FALSE]
    xi0_s    <- xi_0out[keep]
    sig0_s   <- if (model == "gaussian") sigma2_0out[keep] else NULL

    betahat  <- colMeans(beta_s)
    leftci   <- apply(beta_s, 2, stats::quantile, probs = alpha/2)
    rightci  <- apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2)

    res <- list(
      BetaHat = betahat,
      LeftCI = leftci, RightCI = rightci,
      BetaSamples = beta_s,
      LambdaSamples = 1/sqrt(eta_s),
      TauSamples = 1/sqrt(xi0_s),
      Sigma2Samples = sig0_s,
      kappa_0 = kappa_0out[keep]
    )
    return(res)
  }

  ## ---------- MAIN LOOP (TARGET + SOURCES) ----------
  for (i in seq_len(nmc)) {

    if (model == "gaussian") {
      # Gaussian: leave X as-is; residualizes target on current delta
      X_new <- X
      y_new <- c(y0 - drop(X0 %*% delta_new), yA)
      Q     <- Q_gauss
    } else {
      # Logistic: PG-weighted transforms
      X0_w <- sweep(X0, 1, sqrt(omega_0), `*`)
      XA_w <- sweep(XA, 1, sqrt(omega_A), `*`)
      X_new <- rbind(X0_w, XA_w)
      y0_w  <- (y0 - 1/2) / sqrt(omega_0) - drop(X0_w %*% delta_new)
      yA_w  <- (yA - 1/2) / sqrt(omega_A)
      y_new <- c(y0_w, yA_w)
      Q     <- crossprod(X_new)
    }

    # --- Auxiliary/shared coefficients w_k update ---
    wk_samp <- exact_horseshoe_gibbs_step(
      X = X_new, y = y_new, eta = eta_k, xi = xi_k, Q = Q,
      w = w, s = s, p = p, a = a, b = b,
      sigma2 = if (model == "gaussian") sigma2_k else NULL,
      kappa  = kappa_k,
      EB = EB,
      sigma2_0 = if (model == "gaussian") sigma2_0 else 1,  # retained arg for your helper
      N_k = N_k
    )
    wk_new  <- wk_samp$new_beta
    eta_k   <- wk_samp$new_eta
    xi_k    <- wk_samp$new_xi
    if (model == "gaussian") sigma2_k <- wk_samp$new_sigma2

    # --- Target delta update (conditional on wk_new) ---
    if (model == "gaussian") {
      X0_new <- X0
      y0_new <- y0 - drop(X0 %*% wk_new)
      Q0     <- Q0_gauss
    } else {
      # Logistic uses the already weighted design X0_w
      y0_new <- (y0 - 1/2) / sqrt(omega_0) - drop(X0_w %*% wk_new)
      Q0     <- crossprod(X0_w)
      X0_new <- X0_w
    }

    delta_samp <- exact_horseshoe_gibbs_step(
      X = X0_new, y = y0_new, eta = eta_0, xi = xi_0, Q = Q0,
      w = w, s = s, p = p, a = a, b = b,
      sigma2 = if (model == "gaussian") sigma2_0 else NULL,
      kappa  = kappa_0,
      EB = EB
    )
    delta_new <- delta_samp$new_beta
    eta_0     <- delta_samp$new_eta
    xi_0      <- delta_samp$new_xi
    if (model == "gaussian") sigma2_0 <- delta_samp$new_sigma2

    # Compose beta
    beta <- wk_new + delta_new

    # Save draws
    betaout[i, ]  <- beta
    wk_out[i, ]   <- wk_new
    etaout[i, ]   <- eta_0
    eta_kout[i, ] <- eta_k
    xi_0out[i]    <- xi_0
    xi_kout[i]    <- xi_k
    if (model == "gaussian") {
      sigma2_0out[i] <- sigma2_0
      sigma2_kout[i] <- sigma2_k
    }
    kappa_0out[i] <- kappa_0
    kappa_kout[i] <- kappa_k

    # Optional EB step
    if (EB && i == burn + iterEBstep) {
      tmp_k <- eb_em_max(
        kappa = kappa_k, N = N_k,
        xi_out = xi_kout,
        sigma2_out = if (model == "gaussian") sigma2_kout else rep(NA_real_, nmc),
        iterEBstep = iterEBstep, eb_counter = eb_counter, i = i, sir = sir
      )
      kappa_k <- tmp_k$kappa

      tmp_0 <- eb_em_max(
        kappa = kappa_0, N = n0,
        xi_out = xi_0out,
        sigma2_out = if (model == "gaussian") sigma2_0out else rep(NA_real_, nmc),
        iterEBstep = iterEBstep, eb_counter = eb_counter, i = i, sir = sir
      )
      kappa_0 <- tmp_0$kappa

      sigma2_0_eb_est <- if (model == "gaussian") tmp_0$sigma2_est else NULL
      sigma2_k_eb_est <- if (model == "gaussian") tmp_k$sigma2_est else NULL

      eb_counter <- eb_counter + 1
      if (progress_every > 0) {
        val1 <- (exp(kappa_0)/(1+exp(kappa_0))) * p
        val2 <- (exp(kappa_k)/(1+exp(kappa_k))) * p
        message(sprintf("EB updated: kappa_0=%.3f, kappa_k=%.3f, E[sparsity]_0=%.2f, E[sparsity]_k=%.2f",
                        kappa_0, kappa_k, val1, val2))
      }
    }

    # PG updates for logistic
    if (model == "binomial") {
      # sources
      for (j in seq_len(N_k)) {
        omega_A[j] <- pgdraw::pgdraw(1, drop(XA[j, , drop = FALSE] %*% wk_new))
      }
      # target
      for (j in seq_len(n0)) {
        omega_0[j] <- pgdraw::pgdraw(1, drop(X0[j, , drop = FALSE] %*% beta))
      }
    }

    if (progress_every > 0 && i %% progress_every == 0) {
      message(sprintf("Iteration: %d", i))
    }
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
  sig0_s   <- if (model == "gaussian") sigma2_0out[keep] else NULL
  sigk_s   <- if (model == "gaussian") sigma2_kout[keep] else NULL
  kappa0_s <- kappa_0out[keep]
  kappak_s <- kappa_kout[keep]

  betahat         <- colMeans(beta_s)
  betahat_median  <- apply(beta_s, 2, median)
  leftci          <- apply(beta_s, 2, stats::quantile, probs = alpha/2)
  rightci         <- apply(beta_s, 2, stats::quantile, probs = 1 - alpha/2)
  lambda_s        <- 1 / sqrt(eta_s)
  lambda_k_s      <- 1 / sqrt(eta_k_s)
  tau0_s          <- 1 / sqrt(xi0_s)
  tauk_s          <- 1 / sqrt(xik_s)

  res <- list(
    BetaHat = betahat,
    BetaHatMedian = betahat_median,
    LeftCI = leftci, RightCI = rightci,
    # Estimates:
    Sigma2Hat = if (model == "gaussian") mean(sig0_s) else NULL,
    TauHat    = mean(tau0_s),
    LambdaHat = colMeans(lambda_s),
    LambdaKHat = colMeans(lambda_k_s),
    TauKHat   = mean(tauk_s),
    Sigma2KHat = if (model == "gaussian") mean(sigk_s) else NULL,
    # Samples:
    BetaSamples   = beta_s,
    wkSamples     = wk_s,
    LambdaSamples = lambda_s,
    TauSamples    = tau0_s,
    kappa_0 = kappa0_s,
    kappa_k = kappak_s,
    # EB (if computed)
    sigma2_0_eb_est = sigma2_0_eb_est,
    sigma2_k_eb_est = sigma2_k_eb_est,
    model = model
  )

  return(res)
}
