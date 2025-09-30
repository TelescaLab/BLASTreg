#' Empirical Bayes update for global sparsity in Horseshoe via EM / Newton steps
#'
#' @description
#' Maximizes the marginal likelihood of the horseshoe global scale through a
#' short EM-style routine (implemented as safeguarded Newton updates) to produce
#' an updated value of the sparsity hyperparameter `kappa`. Optionally performs
#' a self–importance reweighting (SIR) refinement using the ratio of likelihoods
#' at the updated vs. baseline `kappa`.
#'
#' @details
#' The routine consumes posterior draws of the global shrinkage parameter
#' \eqn{\xi} (and, optionally, \eqn{\sigma^2}) from earlier MCMC iterations and
#' optimizes the marginal log-likelihood with respect to `kappa`, where the
#' working global scale is \eqn{\psi = \exp(\kappa)\sqrt{\sigma^2/N}}.
#'
#' The function calls an internal helper (`emp_bayes_step()`) that:
#'   - computes gradient and Hessian of the objective,
#'   - performs a Newton step with backtracking line search,
#'   - (optionally) applies SIR weights to stabilize/adjust the update.
#'
#' Only the subset of indices corresponding to the most recent EB window is
#' used (from `eb_counter * iterEBstep + 1` through `i`). Convergence is when
#' `|grad| < tol`, or when `maxiter` Newton iterations have elapsed.
#'
#' @param kappa Numeric scalar. Current value of the sparsity hyperparameter to
#'   be updated.
#' @param N Integer. Effective sample size used to scale the global parameter
#'   \eqn{\psi = \exp(\kappa)\sqrt{\sigma^2/N}}.
#' @param xi_out Numeric vector of past/posterior draws of the global parameter
#'   \eqn{\xi}; only the most recent EB window is used internally.
#' @param sigma2_out Numeric vector of past/posterior draws of \eqn{\sigma^2}.
#'   (Note: in the current implementation of `emp_bayes_step()` the EB update
#'   uses a working \eqn{\sigma^2 = 1} for stability; `sigma2_out` is accepted
#'   for future extensions and API symmetry.)
#' @param iterEBstep Integer. Size of the EB window (how many iterations’ worth
#'   of draws to use per EB refresh).
#' @param eb_counter Integer. How many EB refreshes have already occurred; used
#'   to compute the start of the EB window.
#' @param i Integer. Current MCMC iteration (end index of the EB window).
#' @param sir Logical. If `TRUE`, apply a SIR-based refinement using weights
#'   proportional to the likelihood ratio at the updated vs. baseline `kappa`.
#' @param maxiter Integer. Maximum Newton steps per EB call.
#' @param tol Numeric. Convergence tolerance on the absolute gradient.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{kappa}{Updated EB estimate of the sparsity parameter.}
#'   \item{grad}{Final (weighted) gradient at the returned `kappa`.}
#'   \item{it}{Number of Newton iterations performed.}
#' }
#'
#' @seealso
#' Internal helpers: `emp_bayes_step()`, `eb_obj_fun()`, `eb_loglik()`,
#' `eb_grad()`, `eb_hess()`. Higher-level samplers such as `blast_select()`
#' and `blast_oracle()` call this routine to refresh global sparsity.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 200
#' xi_draws <- rexp(1000, rate = 1)   # toy xi draws
#' sig2_draws <- rep(1, 1000)
#' out <- eb_em_max(
#'   kappa = -0.5, N = N,
#'   xi_out = xi_draws, sigma2_out = sig2_draws,
#'   iterEBstep = 200, eb_counter = 0, i = 200,
#'   sir = TRUE
#' )
#' out$kappa
#' }
#'
#' @export
eb_em_max <- function(kappa, N, xi_out, sigma2_out, iterEBstep, eb_counter, i, sir = T, maxiter = 20, tol = 1e-6){
  stop = FALSE
  it = 0
  while(!stop){
    if(it > maxiter) break
    eb <- emp_bayes_step(kappa_old = kappa, N = N, xi_out = xi_out,
                         sigma2_out = sigma2_out,
                         eb_inds = (eb_counter*iterEBstep + 1):(i),
                         compute_grad = T, compute_hess = T, sir = sir)
    stop = abs(eb$grad) < tol
    kappa <- eb$kappa_new
    it = it + 1
  }
  return(list(kappa = kappa, grad = eb$grad, it = it))
}

#' @keywords internal
#' @noRd
emp_bayes_step <- function(kappa_old, N, xi_out, sigma2_out, eb_inds, compute_grad = TRUE, compute_hess = TRUE, sir = FALSE,
                           maxiter = 20, tol = 1e-6){
  # set sigma2 to 1 for empirical bayes step
  sigma2 <- rep(1, length(eb_inds))
  xi <- xi_out[eb_inds]
  alpha = 1
  if(compute_grad){
    grad <- eb_grad(kappa_old, xi, sigma2, N)
  }
  if(compute_hess){
    hess <- eb_hess(kappa_old, xi, sigma2, N)
  } else{
    hess <- 1
  }
  # Newton-Raphson update with backtracking
  kappa_new <- kappa_old - alpha * (grad / hess)
  obj_old <- eb_obj_fun(kappa = kappa_old, N, xi, sigma2)
  obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2)
  while((obj_old - obj_new) > tol){
    alpha = alpha / 2
    kappa_new <- kappa_old - alpha * (grad / hess)
    obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2)
  }
  if(sir){
    i = 0
    lik_kappa_0 <- eb_loglik(kappa_old, N, xi, sigma2, exponentiate = TRUE)
    grad = Inf
    while(grad > tol & i < maxiter){
      if(i > maxiter) break
      alpha = 1
      kappa_k <- kappa_new
      lik_kappa_k <- eb_loglik(kappa_k, N, xi, sigma2, exponentiate = TRUE)
      ws <- lik_kappa_k / lik_kappa_0
      grad <- eb_grad(kappa_k, xi, sigma2, N, ws = ws)
      hess <- eb_hess(kappa_k, xi, sigma2, N, ws = ws)
      obj_old <- eb_obj_fun(kappa = kappa_k, N, xi, sigma2, ws = ws)
      kappa_new <- kappa_k - alpha * (grad / hess)
      obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2, ws = ws)
      j = 0
      while((obj_old - obj_new) > tol){
        alpha = alpha / 2
        kappa_new <- kappa_k - alpha * (grad / hess)
        obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2, ws = ws)
        grad <- eb_grad(kappa_new, xi, sigma2, N, ws = ws)
        if(j > maxiter) break
        j = j + 1
      }
      i = i + 1
    }
  }
  return(list(kappa_new = kappa_new, grad = grad, hess = hess))
}

#' @keywords internal
#' @noRd
eb_obj_fun <- function(kappa, N, xi, sigma2, ws = 1){
  logliks <- eb_loglik(kappa, N, xi, sigma2)
  obj <- mean(logliks * ws)
  return(obj)
}

#' @keywords internal
#' @noRd
eb_loglik <- function(kappa, N, xi, sigma2 = 1, exponentiate = FALSE){
  xi_inv <- 1 / xi
  psi = exp(kappa) * sqrt(sigma2 / N)
  liks <- log(2/pi) - log(psi + 1/(xi * psi))
  if(exponentiate) liks <- exp(liks)
  return(liks)
}

#' @keywords internal
#' @noRd
eb_grad <- function(kappa, xi, sigma2, N, ws = 1){
  xi_inv <- 1 / xi
  u = exp(2*kappa) * (sigma2 / N)
  grads <- -(u - xi_inv) / (u + xi_inv)
  grad <- mean(grads * ws)
}

#' @keywords internal
#' @noRd
eb_hess <- function(kappa, xi, sigma2, N, ws = 1){
  xi_inv <- 1 / xi
  u = exp(2*kappa) * (sigma2 / N)
  hess_vals <- -4*u*xi_inv / (u + xi_inv)^2
  hess <- mean(hess_vals * ws)
}
