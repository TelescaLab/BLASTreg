#' One Gibbs update for the (exact) horseshoe regression
#'
#' @description
#' Performs a single sampling sweep for the **horseshoe** prior in a
#' Gaussian linear model. The step updates:
#' (i) the global shrinkage \eqn{\xi} via a random walk on \eqn{\log \xi},
#' (ii) the regression coefficients \eqn{\beta}, and
#' (iii) the local shrinkage parameters \eqn{\eta_j}.
#'
#' Includes small numerical-stability tweaks (capping extreme \eqn{\eta}/\eqn{\xi}
#' and using \eqn{p<N} vs. \eqn{p\ge N} algebra). Portions adapt code from the
#' **Mhorseshoe** package (Johndrow et al., 2020); see **References**.
#'
#' @details
#' Let \eqn{Q = X^\top X}. When \eqn{p < N}, the \eqn{\beta}-update uses solves with
#' \eqn{Q_\* = \xi\,\mathrm{diag}(\eta) + Q}. When \eqn{p \ge N}, it uses
#' \eqn{M = I_N + X \mathrm{diag}(1/\eta) X^\top / \xi}.
#'
#' The global parameter \eqn{\xi} is proposed on the log scale with variance \code{s}.
#' You can optionally **truncate** the proposal below at \code{log(truncation_threshold)}
#' by setting \code{truncated_global_update = TRUE}. In that case, the MH accept
#' ratio includes the necessary truncated-normal normalizing-constant correction.
#'
#' The global-scale factor \eqn{\psi} is set either from an empirical-Bayes sparsity
#' parameter \code{kappa} when \code{EB = TRUE}, or from a target number of nonzeros
#' \code{p_star} (if supplied) otherwise.
#'
#' @param X Numeric matrix of predictors (\eqn{N \times p}).
#' @param y Numeric response vector (length \eqn{N}).
#' @param eta Numeric vector (length \eqn{p}) of current local scales \eqn{\eta_j}.
#' @param xi Numeric scalar: current global scale \eqn{\xi}.
#' @param Q Precomputed \eqn{X^\top X} (\eqn{p \times p}).
#' @param w Positive scalar: Inverse-Gamma prior shape parameter used in updates.
#' @param s Nonnegative scalar: proposal variance for the random walk on \eqn{\log \xi}.
#' @param a,b Positive scalars: tuning parameters for the rejection sampler of \eqn{\eta}.
#' @param p Integer: number of predictors.
#' @param sigma2 Numeric scalar or `NULL`: error variance. If `NULL`, a default of `1`
#'   is used; otherwise it is sampled from the usual Inverse-Gamma update.
#' @param kappa Numeric: EB sparsity parameter when \code{EB = TRUE}.
#' @param EB Logical: if `TRUE`, use EB scaling via \code{kappa}; otherwise use
#'   \code{p_star} (if provided) or a neutral default.
#' @param p_star Optional integer: target number of nonzeros when \code{EB = FALSE}.
#' @param truncated_global_update Logical: if `TRUE`, propose \eqn{\log \xi} from a
#'   lower-truncated normal with lower bound \code{log(truncation_threshold)}.
#' @param truncation_threshold Numeric lower bound used when
#'   \code{truncated_global_update = TRUE}. Must be > 0.
#'
#' @return A list with:
#' - `new_beta`: numeric vector length \eqn{p} (updated coefficients),
#' - `new_xi`: updated global scale \eqn{\xi},
#' - `new_sigma2`: updated error variance \eqn{\sigma^2},
#' - `new_eta`: vector length \eqn{p} (updated local scales).
#'
#' @note
#' Code structure and several numerical tricks are adapted from
#' \strong{Mhorseshoe} (Johndrow et al., 2020). Extreme values of \eqn{\eta} and
#' \eqn{\xi} are lightly capped for stability in high dimensions.
#'
#' @references
#' Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020).
#' \emph{Mhorseshoe}: Efficient sampling for the horseshoe prior. R package.
#'
#' @seealso
#' \code{\link[Mhorseshoe]{Mhorseshoe}}; higher-level routines such as
#' \code{blast_select()} / \code{blast_oracle()} call this step internally.
#'
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats rnorm runif rgamma pnorm
#' @export
exact_horseshoe_step <- function(
    X, y, eta, xi, Q, w, s, a, b, p, sigma2, kappa, EB, p_star = NULL,
    truncated_global_update = FALSE, truncation_threshold = NULL
) {
  N <- nrow(X)

  # -----------------------------
  # Helpers
  # -----------------------------
  log_trunc_norm_const <- function(mu, lower, sd) {
    # log P(N(mu, sd^2) > lower)
    stats::pnorm((lower - mu) / sd, lower.tail = FALSE, log.p = TRUE)
  }

  # -----------------------------
  # Propose log(xi): optionally truncated RW on log-scale
  # -----------------------------
  if (!is.finite(s) || s < 0) stop("s must be nonnegative (proposal variance on log(xi)).")
  sd_prop <- sqrt(s)

  if (truncated_global_update) {
    if (is.null(truncation_threshold) || !is.finite(truncation_threshold) || truncation_threshold <= 0) {
      stop("truncation_threshold must be finite and > 0 when truncated_global_update = TRUE.")
    }
    lower <- log(truncation_threshold)

    # If s == 0, proposal is degenerate at log(xi); truncation would be incompatible
    # if log(xi) < lower. We handle that cleanly.
    if (s == 0) {
      if (log(xi) < lower) stop("s=0 but current log(xi) is below truncation bound; cannot propose.")
      log_xi_prop <- log(xi)
    } else {
      log_xi_prop <- truncnorm::rtruncnorm(
        1, a = lower, b = Inf, mean = log(xi), sd = sd_prop
      )
    }
  } else {
    log_xi_prop <- if (s == 0) log(xi) else stats::rnorm(1, mean = log(xi), sd = sd_prop)
  }
  xi_prop <- exp(log_xi_prop)

  # -----------------------------
  # Light numerical caps (stability)
  # -----------------------------
  eta_cap <- 1e8
  xi_cap  <- 1e8
  eta <- pmin(eta, eta_cap)
  xi  <- min(xi,  xi_cap)
  xi_prop <- min(xi_prop, xi_cap)

  # -----------------------------
  # Global scale factor psi (EB vs p_star/default)
  # -----------------------------
  if (!EB) {
    if (is.null(p_star)) {
      psi <- 1
    } else {
      if (!is.finite(p_star) || p_star <= 0 || p_star >= p) stop("p_star must satisfy 0 < p_star < p.")
      psi <- (p_star / (p - p_star)) * sqrt(1 / N)
    }
  } else {
    psi <- exp(kappa) * sqrt(1 / N)
  }

  # -----------------------------
  # xi update (Metropolis step on log-scale proposal)
  # -----------------------------
  if (p < N) {
    Xy <- crossprod(X, y)
    y_square <- drop(crossprod(y))

    Q_star <- xi * diag(eta) + Q
    m <- solve(Q_star, Xy)
    ymy <- y_square - drop(crossprod(y, X %*% m))

    if (s != 0) {
      new_Q_star <- xi_prop * diag(eta) + Q
      new_m <- solve(new_Q_star, Xy)
      new_ymy <- y_square - drop(crossprod(y, X %*% new_m))

      # cM corresponds to det(Q_star)/xi^p up to a constant (matches Mhorseshoe trick)
      cM     <- (diag(chol(Q_star))^2)     / xi
      new_cM <- (diag(chol(new_Q_star))^2) / xi_prop

      curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w + ymy) -
        log(sqrt(xi) / psi * (1 + xi * psi^2))
      new_ratio  <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log(sqrt(xi_prop) / psi * (1 + xi_prop * psi^2))

      log_accept <- (new_ratio - curr_ratio) + (log(xi_prop) - log(xi))  # Jacobian for xi = exp(log_xi)

      # Truncated RW correction: proposal q(log_xi_prop | log_xi) not symmetric
      if (truncated_global_update) {
        lower <- log(truncation_threshold)
        logZ_curr <- log_trunc_norm_const(mu = log(xi),     lower = lower, sd = sd_prop)
        logZ_new  <- log_trunc_norm_const(mu = log(xi_prop), lower = lower, sd = sd_prop)
        log_accept <- log_accept + (logZ_curr - logZ_new)
      }

      if (log(stats::runif(1)) < min(0, log_accept)) {
        xi    <- xi_prop
        ymy   <- new_ymy
        Q_star <- new_Q_star
      }
    }

  } else {
    DX  <- (1/eta) * t(X)
    XDX <- X %*% DX
    M <- diag(N) + XDX/xi
    m <- solve(M, y)
    ymy <- drop(crossprod(y, m))

    if (s != 0) {
      new_M <- diag(N) + XDX/xi_prop
      new_m <- solve(new_M, y)
      new_ymy <- drop(crossprod(y, new_m))

      cM     <- diag(chol(M))^2
      new_cM <- diag(chol(new_M))^2

      curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w + ymy) -
        log((sqrt(xi) / psi) * (1 + xi * psi^2))
      new_ratio  <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log((sqrt(xi_prop) / psi) * (1 + xi_prop * psi^2))

      log_accept <- (new_ratio - curr_ratio) + (log(xi_prop) - log(xi))

      if (truncated_global_update) {
        lower <- log(truncation_threshold)
        logZ_curr <- log_trunc_norm_const(mu = log(xi),     lower = lower, sd = sd_prop)
        logZ_new  <- log_trunc_norm_const(mu = log(xi_prop), lower = lower, sd = sd_prop)
        log_accept <- log_accept + (logZ_curr - logZ_new)
      }

      if (log(stats::runif(1)) < min(0, log_accept)) {
        xi  <- xi_prop
        ymy <- new_ymy
        M   <- new_M
      }
    }
  }

  # -----------------------------
  # sigma^2 update
  # -----------------------------
  if (is.null(sigma2)) {
    sigma2 <- 1
  } else {
    sigma2 <- 1 / stats::rgamma(1, shape = (w + N)/2, rate = (w + ymy)/2)
  }

  # -----------------------------
  # beta update
  # -----------------------------
  diag_D <- 1 / (eta * xi)
  u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
  f <- stats::rnorm(n = N, mean = 0, sd = 1)
  v <- X %*% u + f
  U <- diag_D * t(X)

  if (p < N) {
    yv   <- y / sqrt(sigma2) - v
    Xyv  <- crossprod(X, yv)
    m    <- solve(Q_star, Xyv)
    m_st <- yv - X %*% m
    new_beta <- sqrt(sigma2) * (u + U %*% m_st)
  } else {
    v_st <- solve(M, (y / sqrt(sigma2) - v))
    new_beta <- sqrt(sigma2) * (u + U %*% v_st)
  }

  # -----------------------------
  # eta update (local scales)
  # -----------------------------
  eta <- Mhorseshoe:::rejection_sampler((new_beta^2) * xi / (2 * sigma2), a, b)
  eta <- pmax(eta, .Machine$double.eps)

  list(new_beta = new_beta, new_xi = xi, new_sigma2 = sigma2, new_eta = eta)
}
