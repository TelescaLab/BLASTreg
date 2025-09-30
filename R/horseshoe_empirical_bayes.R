#' Empirical Bayes maximization via EM algorithm
#'
#' @param kappa initial value of global shrinkage prior parameter
#'
#' @param N total number of observations
#' @param xi_out global shrinkage parameter draws from posterior
#' @param sigma2_out error variance draws from posterior
#' @param iterEBstep number of iterations for empirical bayes step
#' @param eb_counter counter for empirical bayes step
#' @param i current iteration
#' @param sir whether to use SIR for empirical bayes
#' @param maxiter maximum number of iterations for EB step
#' @param tol tolerance for EB step
#'
#' @export
eb_em_max <- function(kappa, N, xi_out, sigma2_out, iterEBstep, eb_counter, i, sir = F, maxiter = 20, tol = 1e-6){
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


emp_bayes_step <- function(kappa_old, N, xi_out, sigma2_out, eb_inds, compute_grad = T, compute_hess = T, sir = F,
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
  # Newton-Raphson update
  kappa_new <- kappa_old - alpha * (grad / hess)
  # Ensure step size is appropriate and is increasing the objective
  obj_old <- eb_obj_fun(kappa = kappa_old, N, xi, sigma2)
  obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2)
  while((obj_old - obj_new) > tol){
    alpha = alpha / 2
    kappa_new <- kappa_old - alpha * (grad / hess)
    obj_new <- eb_obj_fun(kappa = kappa_new, N, xi, sigma2)
  }
 if(sir){
    i = 0
    lik_kappa_0 <- eb_loglik(kappa_old, N, xi, sigma2, exponentiate = T)
    grad = Inf
    while(grad > tol & i < maxiter){
      if(i > maxiter) break
      alpha = 1
      kappa_k <- kappa_new
      lik_kappa_k <- eb_loglik(kappa_k, N, xi, sigma2, exponentiate = T)
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

eb_obj_fun <- function(kappa, N, xi, sigma2, ws = 1){
  logliks <- eb_loglik(kappa, N, xi, sigma2)
  obj <- mean(logliks * ws)
  return(obj)
}

eb_loglik <- function(kappa, N, xi, sigma2 = 1, exponentiate = F){
  xi_inv <- 1 / xi
  psi = exp(kappa) * sqrt(sigma2 / N)
  liks <- log(2/pi) - log(psi + 1/(xi * psi))
  if(exponentiate) liks <- exp(liks)
  return(liks)
}

eb_grad <- function(kappa, xi, sigma2, N, ws = 1){
  xi_inv <- 1 / xi
  u = exp(2*kappa) * (sigma2 / N)
  grads <- -(u - xi_inv) / (u + xi_inv)
  grad <- mean(grads * ws)
}

eb_hess <- function(kappa, xi, sigma2, N, ws = 1){
  xi_inv <- 1 / xi
  u = exp(2*kappa) * (sigma2 / N)
  hess_vals <- -4*u*xi_inv / (u + xi_inv)^2
  hess <- mean(hess_vals * ws)
}


