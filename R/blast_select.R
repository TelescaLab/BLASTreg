#' Bayesian Transfer Learning Shrinkage Regression with Unknown Informative Samples
#'
#' @description
#' Runs an MCMC sampler for a BLAST linear regression with an unknown set
#' of informative source studies. The algorithm alternates between updating
#' (i) target bias coefficients, (ii) auxiliary coefficients for the currently
#' informative and non-informative sets, and (iii) the inclusion vector
#' \eqn{\gamma} via a marginal-likelihood–based update.
#'
#' @details
#' The design matrix `X` and response `y` should be stacked as the target study
#' followed by the auxiliary studies according to `n.vec`. The first element of
#' `n.vec` is the target sample size; the remaining `K = length(n.vec) - 1`
#' entries correspond to auxiliary studies in sequence.
#'
#' @param X Predictor matrix of dimension \eqn{(\sum n_k) \times p}.
#' @param y Response vector of length \eqn{\sum n_k}.
#' @param n.vec Integer vector of study sample sizes `c(n0, n1, ..., nK)`.
#'   The first entry `n0` is the target size; the rest are auxiliary-study sizes.
#' @param burn Number of burn-in iterations.
#' @param iter Number of post–burn-in iterations to keep.
#' @param a Parameter for the rejection sampler used in horseshoe updates.
#' @param b Parameter for the rejection sampler used in horseshoe updates.
#' @param s Tuning parameter for the proposal distribution in the horseshoe step.
#' @param tau Global shrinkage parameter (prior scale; used for initialization).
#' @param sigma2 Error variance (used for initialization).
#' @param w Inverse-Gamma prior shape parameter used in horseshoe updates.
#' @param alpha Credible-interval level (e.g., `0.05` for 95% CIs).
#' @param iterEBstep Number of empirical Bayes EM steps to estimate global shrinkage scale parameter; if `0`, empirical Bayes is disabled.
#' @param gamma_init Optional binary vector of length `K` giving the initial
#'   informative-set indicator; if `NULL`, a candidate set is chosen via
#'   `construct_candidate_set()`.
#' @param temp_scale Numeric temperature multiplier for \eqn{\gamma}-update
#'   probabilities (values > 1 flatten differences; < 1 sharpen). If `NULL`,
#'   the function initializes `temp_scale <- 1/p` and enables adaptive tempering.
#' @param enforce_delta_stronger_shrinkage Logical; if TRUE, the global proposal
#'   for the target “contrast” block (delta) is lower-truncated at the current
#'   global scale of the informative/source block, enforcing xi_delta ≥ xi_source.
#'
#' @return
#' A list with posterior summaries and draws:
#'
#' - `BetaHat`, `BetaHatMedian`: posterior mean/median of \eqn{\beta}.
#' - `LeftCI`, `RightCI`: marginal \(100(1-\alpha)\%\) credible bounds.
#' - `Sigma2Hat`, `TauHat`, `LambdaHat`: posterior means of \eqn{\sigma^2},
#'   \eqn{\tau}, and local shrinkage \eqn{\lambda_j}.
#' - `BetaSamples`, `LambdaSamples`, `TauSamples`, `Sigma2Samples`:
#'   post–burn-in MCMC draws.
#' - `GammaSamples`: draws of the informative-set indicator \eqn{\gamma}.
#' - `W0Samples`, `WA_Samples`: auxiliary coefficients (non-informative vs informative).
#' - `GammaSums`: inclusion counts for each auxiliary study across kept draws.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- simulate_multistudy_regression(
#'   p = 50, s = 3, M = 5, size.A0 = 2,
#'   n.vec = c(100, rep(100, 5)),
#'   sig.beta = 0.5, sig.delta1 = 0.3, sig.delta2 = 1,
#'   contam_pct = 0.01, type = "gaussian"
#' )
#'
#' fit <- blast_select(
#'   X = df$X, y = df$y, n.vec = df$n.vec,
#'   burn = 100, iter = 200, adapt_temp_scale = TRUE
#' )
#' str(fit$BetaHat)
#' }
#'
#' @export
blast_select <- function(
    X, y, n.vec,
    burn = 1000, iter = 3000,
    a = 1/5, b = 10,
    s = 0.8, tau = 1, sigma2 = 1, w = 1, alpha = 0.05,
    iterEBstep = 0,
    gamma_init = NULL,
    temp_scale = NULL,
    enforce_delta_stronger_shrinkage = FALSE
) {
  ## ----- Basic dimensions & bookkeeping -----
  p <- ncol(X)
  N <- sum(n.vec)
  K <- length(n.vec) - 1

  ## Adaptive-tempering counters (used only if adapt_temp_scale = TRUE)
  accepted_count <- 0L
  total_gamma_proposals <- 0L
  target_acceptance_rate <- 0.3

  ## ----- Initialize gamma (informative-set indicator) -----
  if (is.null(gamma_init)) {
    candidate_sets <- construct_candidate_set(X = X, y = y, n.vec = n.vec)
    gamma <- rep(0, K)
    gamma[candidate_sets[[1]]] <- 1
  } else {
    gamma <- gamma_init
  }

  if(is.null(temp_scale)){
    temp_scale <- 1 / p
    adapt_temp_scale <- TRUE
  } else{
    adapt_temp_scale <- FALSE
  }

  ## ----- Partition target vs. auxiliary -----
  n0 <- n.vec[1]
  N_k <- N - n0
  X_0 <- X[1:n0, , drop = FALSE]
  y_0 <- y[1:n0]

  ## Target prior state
  eta_0 <- rep(1, p)
  xi_0  <- tau^(-2)
  Q0    <- crossprod(X_0)
  sigma2_0 <- 1
  kappa_0  <- -0.5

  ## Auxiliary prior state
  etaA      <- rep(1, p)
  xiA       <- tau^(-2)
  eta_A_bar <- rep(1, p)
  xi_A_bar  <- tau^(-2)
  wA_new    <- rep(0, p)

  ## Gamma prior
  gamma_prior_probs <- rep(0.5, K)

  ## Auxiliary likelihood variances & kappas
  sigma2_A      <- 1
  sigma2_A_bar  <- 1
  kappa_A       <- -0.5
  kappa_A_bar   <- -0.5

  ## Storage
  nmc <- burn + iter
  xi_0out          <- rep(0, nmc)
  xiA_out          <- rep(0, nmc)
  xiA_bar_out      <- rep(0, nmc)
  sigma2_0out      <- rep(0, nmc)
  sigma2_Aout      <- rep(0, nmc)
  sigma2_A_bar_out <- rep(0, nmc)
  betaout          <- matrix(0, nrow = nmc, ncol = p)
  wAout            <- matrix(0, nrow = nmc, ncol = p)
  w0out            <- matrix(0, nrow = nmc, ncol = p)
  deltaout         <- matrix(0, nrow = nmc, ncol = p)
  gammaout         <- matrix(0, nrow = nmc, ncol = K)
  etaout           <- matrix(0, nrow = nmc, ncol = p)

  ## Helpful aliases for target sufficient stats
  X0_X0_t <- Q0
  X0_y0   <- crossprod(X_0, y_0)

  ## Precompute per-study XtX and Xty for auxiliaries
  precomp <- precompute_XtXXty(X, y, K, n.vec)

  ## ----- Optional empirical Bayes (target-only warm-start) -----
  EB <- iterEBstep != 0
  if (EB) {
    EB_burn <- 1000
    EB_step_total <- iterEBstep + EB_burn

    eb_beta_samples   <- matrix(0, nrow = EB_step_total + 1, ncol = p)
    eb_eta_samples    <- matrix(0, nrow = EB_step_total + 1, ncol = p)
    eb_xi_samples     <- rep(0, EB_step_total + 1)
    eb_sigma2_samples <- rep(0, EB_step_total + 1)

    eb_beta_samples[1, ] <- 0
    eb_eta_samples[1, ]  <- 1
    eb_xi_samples[1]     <- xi_0
    eb_sigma2_samples[1] <- sigma2_0

    for (j in 2:(EB_step_total + 1)) {
      targ_sample <- exact_horseshoe_gibbs_step(
        X_0, y_0,
        eta = eb_eta_samples[j - 1, ], xi = eb_xi_samples[j - 1],
        Q = Q0, w = w, s = s, a = a, b = b, p = p,
        sigma2 = eb_sigma2_samples[j - 1],
        kappa = kappa_0, EB = FALSE
      )
      eb_beta_samples[j, ]   <- targ_sample$new_beta
      eb_eta_samples[j, ]    <- targ_sample$new_eta
      eb_xi_samples[j]       <- targ_sample$new_xi
      eb_sigma2_samples[j]   <- targ_sample$new_sigma2
    }

    tmp <- eb_em_max(
      kappa = kappa_0, N = n0,
      xi_out = eb_xi_samples[(EB_burn + 1):(EB_step_total + 1)],
      sigma2_out = eb_sigma2_samples[(EB_burn + 1):(EB_step_total + 1)],
      iterEBstep = iterEBstep, eb_counter = 0,
      i = iterEBstep
    )

    kappa_A     <- tmp$kappa
    kappa_0     <- tmp$kappa
    kappa_A_bar <- tmp$kappa

  }
  ## =======================
  ## =       MCMC         =
  ## =======================
  for (i in 1:nmc) {

    ## (1) Set indices based on current gamma
    data_inds <- set_data_inds(gamma, n.vec)
    ind.kA <- data_inds$ind.kA
    ind.ni <- data_inds$ind.ni

    ## (2) Update target delta given current wA_new
    delta_samp <- exact_horseshoe_gibbs_step(
      X = X_0,
      y = y_0 - X_0 %*% wA_new,
      eta = eta_0, xi = xi_0, Q = Q0,
      w = w, s = s, a = a, b = b, p = p,
      sigma2 = sigma2_0, kappa = kappa_0, EB = EB,
      truncated_global_update = enforce_delta_stronger_shrinkage,
      truncation_threshold    = if (enforce_delta_stronger_shrinkage) xiA else NULL
    )

    delta_new <- delta_samp$new_beta
    xi_0      <- delta_samp$new_xi
    sigma2_0  <- delta_samp$new_sigma2
    eta_0     <- delta_samp$new_eta

    ## (3) Update auxiliary coefficients for informative set
    if (data_inds$informative) {
      X_A  <- X[ind.kA, , drop = FALSE]
      y_A  <- y[ind.kA]
      newX <- rbind(X_0, X_A)
      newy <- c(y_0 - X_0 %*% delta_new, y_A)
      Q    <- crossprod(newX)
    } else {
      newX <- X_0
      newy <- y_0 - X_0 %*% delta_new
      Q    <- crossprod(newX)
    }

    wA_samp <- exact_horseshoe_gibbs_step(
      X = newX, y = newy,
      eta = etaA, xi = xiA, Q = Q,
      w = w, s = s, a = a, b = b, p = p,
      sigma2 = sigma2_A, kappa = kappa_A, EB = EB,
      truncated_global_update = FALSE,
      truncation_threshold = NULL
    )

    wA_new   <- wA_samp$new_beta
    xiA      <- wA_samp$new_xi
    sigma2_A <- wA_samp$new_sigma2
    etaA     <- wA_samp$new_eta
    wAout[i, ] <- wA_new

    ## (4) Update auxiliary coefficients for non-informative set
    if (data_inds$noninformative) {
      X_A_bar <- X[ind.ni, , drop = FALSE]
      y_A_bar <- y[ind.ni]
      Q_A_bar <- crossprod(X_A_bar)
    } else {
      pseudo_inds <- sample((n0 + 1):N, size = floor(0.03 * N_k))
      y_A_bar     <- y[sample(pseudo_inds)]
      X_A_bar     <- X[pseudo_inds, , drop = FALSE]
      Q_A_bar     <- crossprod(X_A_bar)
    }

    w0_samp <- exact_horseshoe_gibbs_step(
      X = X_A_bar, y = y_A_bar,
      eta = eta_A_bar, xi = xi_A_bar, Q = Q_A_bar,
      w = w, s = s, a = a, b = b, p = p,
      sigma2 = sigma2_A_bar, kappa = kappa_A_bar, EB = EB,
      truncated_global_update = FALSE,
      truncation_threshold = NULL
    )

    w0_new       <- w0_samp$new_beta
    xi_A_bar     <- w0_samp$new_xi
    sigma2_A_bar <- w0_samp$new_sigma2
    eta_A_bar    <- w0_samp$new_eta
    w0out[i, ]   <- w0_new

    ## (5) Compose beta and save partial draws
    deltaout[i, ] <- delta_new
    betaout[i, ]  <- wA_new + delta_new

    ## (6) Diagonal prior terms for marginal-likelihood γ update
    Sigma_A     <- (1 / xiA)      * diag(1 / etaA)
    Sigma_delta <- (1 / xi_0)     * diag(1 / eta_0)
    Sigma_A_bar <- (1 / xi_A_bar) * diag(1 / eta_A_bar)

    ## (7) Update gamma one coordinate at a time
    for (k in 1:K) {
      gamma_tmp <- gamma

      ## Proposal: gamma_k = 1
      gamma_tmp[k] <- 1
      inds.tmp <- set_data_inds(gamma_tmp, n.vec)
      if (!inds.tmp$noninformative) {
        yAb <- 0
        XAb <- t(rep(0, p))
      } else {
        yAb <- y[inds.tmp$ind.ni]
        XAb <- X[inds.tmp$ind.ni, , drop = FALSE]
      }

      aux_num_ind.kA <- as.logical(gamma_tmp)
      aux_num_ind.ni <- !aux_num_ind.kA
      precomp_XtX_A     <- precomp$list_XtX[aux_num_ind.kA]
      precomp_Xty_A     <- precomp$list_Xty[aux_num_ind.kA]
      precomp_XtX_Ab    <- precomp$list_XtX[aux_num_ind.ni]
      precomp_Xty_Ab    <- precomp$list_Xty[aux_num_ind.ni]

      log_p1 <- log_marginal_likelihood(
        y0 = y_0, yA = y[inds.tmp$ind.kA], yAb = yAb,
        X0 = X_0, XA = X[inds.tmp$ind.kA, , drop = FALSE], XAb = XAb,
        d_A = diag(Sigma_A), d_Ab = diag(Sigma_A_bar), d_delta = diag(Sigma_delta),
        X0_X0_t = X0_X0_t, X0_y0 = X0_y0,
        precomp_XtX_A = precomp_XtX_A, precomp_Xty_A = precomp_Xty_A,
        precomp_XtX_Ab = precomp_XtX_Ab, precomp_Xty_Ab = precomp_Xty_Ab,
        informative = inds.tmp$informative, noninformative = inds.tmp$noninformative
      )

      ## Proposal: gamma_k = 0
      gamma_tmp[k] <- 0
      aux_num_ind.kA <- as.logical(gamma_tmp)
      aux_num_ind.ni <- !aux_num_ind.kA
      precomp_XtX_A     <- precomp$list_XtX[aux_num_ind.kA]
      precomp_Xty_A     <- precomp$list_Xty[aux_num_ind.kA]
      precomp_XtX_Ab    <- precomp$list_XtX[aux_num_ind.ni]
      precomp_Xty_Ab    <- precomp$list_Xty[aux_num_ind.ni]
      inds.tmp <- set_data_inds(gamma_tmp, n.vec)

      if (!inds.tmp$informative) {
        yA_prop <- 0
        XA_prop <- t(rep(0, p))
      } else {
        yA_prop <- y[inds.tmp$ind.kA]
        XA_prop <- X[inds.tmp$ind.kA, , drop = FALSE]
      }

      log_p0 <- log_marginal_likelihood(
        y0 = y_0, yA = yA_prop, yAb = y[inds.tmp$ind.ni],
        X0 = X_0, XA = XA_prop, XAb = X[inds.tmp$ind.ni, , drop = FALSE],
        d_A = diag(Sigma_A), d_Ab = diag(Sigma_A_bar), d_delta = diag(Sigma_delta),
        X0_X0_t = X0_X0_t, X0_y0 = X0_y0,
        precomp_XtX_A = precomp_XtX_A, precomp_Xty_A = precomp_Xty_A,
        precomp_XtX_Ab = precomp_XtX_Ab, precomp_Xty_Ab = precomp_Xty_Ab,
        informative = inds.tmp$informative, noninformative = inds.tmp$noninformative
      )

      ## Convert to probabilities (with current temp_scale) & sample gamma_k
      probs <- normalize_log_probabilities(log_p1 = log_p1, log_p0 = log_p0, scale = temp_scale)
      p1 <- probs[1]
      gamma_old <- gamma[k]
      gamma[k] <- rbinom(1, 1, p1)

      ## Track acceptance for adaptive tempering
      if (adapt_temp_scale) {
        total_gamma_proposals <- total_gamma_proposals + 1L
        if (gamma[k] != gamma_old) accepted_count <- accepted_count + 1L
      }
    }

    ## Save iteration-level outputs
    gammaout[i, ]        <- gamma
    etaout[i, ]          <- eta_0
    xi_0out[i]           <- xi_0
    xiA_out[i]           <- xiA
    xiA_bar_out[i]       <- xi_A_bar
    sigma2_0out[i]       <- sigma2_0
    sigma2_Aout[i]       <- sigma2_A
    sigma2_A_bar_out[i]  <- sigma2_A_bar

    ## Adaptive tempering update (small, bounded, diminishing step)
    if (adapt_temp_scale && total_gamma_proposals > 0L) {
      acceptance_prob <- accepted_count / total_gamma_proposals
      step_size <- 0.1 / sqrt(i)                # diminishing step
      temp_scale <- temp_scale + step_size * (acceptance_prob - target_acceptance_rate)
      temp_scale <- min(max(temp_scale, 1 / p^2), 10)  # clamp
      # lightweight progress (comment out if too chatty)
      # cat("Temp scale:", round(temp_scale, 4), " | acc:", round(acceptance_prob, 3), "\n")
    }

  }

  ## ----- Post-processing: drop burn-in & summarize -----
  kept_idx     <- (burn + 1):nmc
  betaout      <- betaout[kept_idx, , drop = FALSE]
  gammaout     <- gammaout[kept_idx, , drop = FALSE]
  lambdaout    <- 1 / sqrt(etaout[kept_idx, , drop = FALSE])
  tauout       <- 1 / sqrt(xi_0out[kept_idx])
  sigma2_0keep <- sigma2_0out[kept_idx]

  betahat         <- colMeans(betaout)
  betahat_median  <- apply(betaout, 2, median)
  lambdahat       <- colMeans(lambdaout)
  tauhat          <- mean(tauout)
  sigma2hat       <- mean(sigma2_0keep)
  gammasums       <- colSums(gammaout)
  leftci          <- apply(betaout, 2, stats::quantile, probs = alpha / 2)
  rightci         <- apply(betaout, 2, stats::quantile, probs = 1 - alpha / 2)

  result <- list(
    BetaHat         = betahat,
    BetaHatMedian   = betahat_median,
    LeftCI          = leftci,
    RightCI         = rightci,
    Sigma2Hat       = sigma2hat,
    TauHat          = tauhat,
    LambdaHat       = lambdahat,
    BetaSamples     = betaout,
    LambdaSamples   = lambdaout,
    TauSamples      = tauout,
    Sigma2Samples   = sigma2_0keep,
    GammaSamples    = gammaout,
    W0Samples       = w0out,
    WA_Samples      = wAout,
    GammaSums       = gammasums
  )

  return(result)
}
