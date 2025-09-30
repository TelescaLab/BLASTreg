#' One Gibbs update for the (exact) horseshoe regression
#'
#' @description
#' Performs a single Gibbs-sampling sweep for the **horseshoe** prior in a
#' Gaussian (or latent-Gaussian) linear model. The step updates:
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
#' by setting \code{truncated_global_update = TRUE}.
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
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 200; p <- 50
#' X <- matrix(rnorm(N * p), N, p)
#' beta <- c(rnorm(5, 0, 2), rep(0, p - 5))
#' y <- as.numeric(X %*% beta + rnorm(N))
#'
#' eta <- rep(1, p); xi <- 1; w <- 1; s <- 0.5; a <- 0.5; b <- 0.5
#' Q <- crossprod(X)
#'
#' out <- exact_horseshoe_gibbs_step(
#'   X, y, eta = eta, xi = xi, Q = Q, w = w, s = s,
#'   a = a, b = b, p = p, sigma2 = 1, kappa = 0, EB = FALSE
#' )
#' str(out)
#' }
#'
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats rnorm runif rgamma
#' @export
exact_horseshoe_gibbs_step <- function(
    X, y, eta, xi, Q, w, s, a, b, p, sigma2, kappa, EB, p_star = NULL,
    truncated_global_update = FALSE, truncation_threshold = NULL
) {
  N <- nrow(X)

  # -- propose log xi (optionally truncated) --
  if (truncated_global_update) {
    if (is.null(truncation_threshold) || truncation_threshold <= 0)
      stop("truncation_threshold must be > 0 when truncated_global_update = TRUE.")
    log_xi <- truncnorm::rtruncnorm(
      1, a = log(truncation_threshold), b = Inf, mean = log(xi), sd = sqrt(s)
    )
  } else {
    log_xi <- stats::rnorm(1, mean = log(xi), sd = sqrt(s))
  }
  new_xi <- exp(log_xi)

  # light numerical caps (stability)
  eta_cap    <- 1e6
  xi_cap     <- 1e6
  new_xi_cap <- 1e6
  eta    <- pmin(eta,    eta_cap)
  xi     <- pmin(xi,     xi_cap)
  new_xi <- pmin(new_xi, new_xi_cap)

  # -- global scale factor psi (EB vs p_star/default) --
  if (!EB) {
    if (is.null(p_star)) {
      psi <- 1
    } else {
      psi <- (p_star / (p - p_star)) * sqrt(1 / N)
    }
  } else {
    psi <- exp(kappa) * sqrt(1 / N)
  }

  # -- xi update (Metropolis on log-scale proposal) --
  if (p < N) {
    Xy <- t(X) %*% y
    y_square <- t(y) %*% y
    Q_star <- xi * diag(eta) + Q
    m <- solve(Q_star, Xy)
    ymy <- y_square - t(y) %*% X %*% m

    if (s != 0) {
      new_Q_star <- new_xi * diag(eta) + Q
      new_m <- solve(new_Q_star, Xy)
      new_ymy <- y_square - t(y) %*% X %*% new_m

      cM     <- (diag(chol(Q_star))^2)     / xi
      new_cM <- (diag(chol(new_Q_star))^2) / new_xi

      curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w + ymy) -
        log(sqrt(xi)     / psi * (1 + xi     * psi^2))
      new_ratio  <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log(sqrt(new_xi) / psi * (1 + new_xi * psi^2))

      acc <- exp(new_ratio - curr_ratio + log(new_xi) - log(xi))
      if (stats::runif(1) < acc) {
        xi    <- new_xi
        ymy   <- new_ymy
        Q_star <- new_Q_star
      }
    }

  } else {
    DX  <- (1/eta) * t(X)
    XDX <- X %*% DX
    M <- diag(N) + XDX/xi
    m <- solve(M, y)
    ymy <- t(y) %*% m

    if (s != 0) {
      new_M <- diag(N) + XDX/new_xi
      new_m <- solve(new_M, y)
      new_ymy <- t(y) %*% new_m

      cM     <- diag(chol(M))^2
      new_cM <- diag(chol(new_M))^2

      curr_ratio <- -sum(log(cM))/2    - ((N + w)/2) * log(w + ymy) -
        log((sqrt(xi)     / psi) * (1 + xi     * psi^2))
      new_ratio  <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log((sqrt(new_xi) / psi) * (1 + new_xi * psi^2))

      acc <- exp(new_ratio - curr_ratio + log(new_xi) - log(xi))
      if (stats::runif(1) < acc) {
        xi  <- new_xi
        ymy <- new_ymy
        M   <- new_M
      }
    }
  }

  # -- sigma^2 update --
  if (is.null(sigma2)) {
    sigma2 <- 1
  } else {
    sigma2 <- 1 / stats::rgamma(1, shape = (w + N)/2, rate = (w + ymy)/2)
  }

  # -- beta update --
  diag_D <- 1 / (eta * xi)
  u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
  f <- stats::rnorm(n = N, mean = 0, sd = 1)
  v <- X %*% u + f
  U <- diag_D * t(X)

  if (p < N) {
    yv   <- (y) / sqrt(sigma2) - v
    Xyv  <- t(X) %*% yv
    m    <- solve(Q_star, Xyv)
    m_st <- yv - X %*% m
    new_beta <- sqrt(sigma2) * (u + U %*% m_st)
  } else {
    v_st <- solve(M, (y / sqrt(sigma2) - v))
    new_beta <- sqrt(sigma2) * (u + U %*% v_st)
  }

  # -- eta update (local scales) --
  eta <- Mhorseshoe:::rejection_sampler((new_beta^2) * xi / (2 * sigma2), a, b)
  eta <- ifelse(eta <= .Machine$double.eps, .Machine$double.eps, eta)

  list(new_beta = new_beta, new_xi = xi, new_sigma2 = sigma2, new_eta = eta)
}
