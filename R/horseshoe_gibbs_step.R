#' One step of the Gibbs sampler for the horseshoe model
#'
#' @param X Predictor matrix
#'
#' @param y Response vector
#' @param eta current value of local shrinkage parameter
#' @param xi current value of global shrinkage parameter
#' @param Q XtX
#' @param w inverse gamma prior shape parameter
#' @param s tuning parameter for the proposal distribution
#' @param a parameter for the rejection sampler
#' @param b parameter for the rejection sampler
#' @param p number of predictors
#' @param sigma2 error variance
#' @param kappa parameter for global shrinkage prior
#' @param EB logical, whether to use empirical Bayes or not
#' @param p_star number of non-zero coefficients if want to specify
#' @param kappa_0_trans logical for whether to transform overall sparsity of bias parameters for TL
#' @param kappa_A_bar_trans logical for whether to transform overall sparsity of noninformative studies for TL
#' @param c0 factor for number of signals for bias
#' @param c_bar factor for number of signals for noninformative studies
#'
#' @export
exact_horseshoe_gibbs_step <- function(X, y, eta, xi, Q, w, s, a, b, p, sigma2, sigma2_0 = NULL, kappa, n0 = NULL, N_k = NULL, EB, p_star = NULL,
                                       kappa_0_trans = FALSE, kappa_A_bar_trans = FALSE, c0 = 1/2, c_bar = 1,
                                       xiA_trunc = FALSE, xiA_trunc_val = NULL){
  N <- nrow(X)
  if(xiA_trunc){
    # Propose log_xi from truncated normal (log(xiA_trunc_val), Inf)
    log_xi <- truncnorm::rtruncnorm(1, a = log(xiA_trunc_val), b = Inf, mean = log(xi), sd = sqrt(s))
  } else{
    log_xi <- stats::rnorm(1, mean = log(xi), sd = sqrt(s))
    #log_xi <- truncnorm::rtruncnorm(1, a = log(1), b = log(p^2), mean = log(xi), sd = sqrt(s))
  }
  # log_xi <- truncnorm::rtruncnorm(1, a = 0, b = log(sqrt(p)), mean = log(xi), sd = sqrt(s))

  new_xi <- exp(log_xi)
  eta_cap <- 1e6  # Example threshold, adjust as needed
  eta <- pmin(eta, eta_cap)
  xi_cap <- 1e6  # Example threshold, adjust as needed
  xi <- pmin(xi, xi_cap)
  new_xi_cap <- 1e6  # Example threshold, adjust as needed
  new_xi <- pmin(new_xi, new_xi_cap)
  if(!EB){
    if(is.null(p_star)){
      psi <- 1
    } else{
      psi <- (p_star / (p - p_star)) * sqrt(1 / N)
    }

  } else{
    psi <- exp(kappa) * sqrt(1 / N)
    # if(kappa_0_trans){
    #   # c is the factor for
    #   psi <- c0 * exp(kappa) / (1 + exp(kappa)*(1 - c0)) * sqrt(sigma2 / N)
    # } else if(kappa_A_bar_trans){
    #   psi <- c_bar * exp(kappa) / (1 + exp(kappa)*(1 - c_bar)) * sqrt(sigma2 / N)
    # } else{
    #   psi <- exp(kappa) * sqrt(sigma2 / N)
    # }
  }
  # xi update
  if (p < N) {
    Xy <- t(X) %*% y
    y_square <- t(y) %*% y
    Q_star <- xi * diag(eta) + Q
    m <- tryCatch(
      solve(Q_star, Xy),
      error = function(e) {
        message("Singularity detected at iteration browser() will launch...")
        browser()   # <-- will drop you into interactive debugging
      }
    )
    #m <- solve(Q_star, Xy)
    ymy <- y_square - t(y) %*% X %*% m
    if (s != 0) {
      new_Q_star <- new_xi * diag(eta) + Q
      #new_m <- solve(new_Q_star, Xy)
      new_m <- tryCatch(
        solve(new_Q_star, Xy),
        error = function(e) {
          message("Singularity detected at iteration browser() will launch...")
          browser()   # <-- will drop you into interactive debugging
        }
      )
      new_ymy <- y_square - t(y) %*% X %*% new_m
      cM <- (diag(chol(Q_star))^2)/xi
      new_cM <- (diag(chol(new_Q_star))^2)/new_xi
      curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w+ymy) -
        log(sqrt(xi) / psi * (1 + xi * psi^2))
      new_ratio <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log(sqrt(new_xi) / psi * (1 + new_xi * psi^2))
      acceptance_probability <- exp(new_ratio - curr_ratio + log(new_xi) -
                                      log(xi))
      u <- stats::runif(n = 1, min = 0, max = 1)
      if (u < acceptance_probability) {
        xi <- new_xi
        ymy <- new_ymy
        Q_star <- new_Q_star
      }
    }
  } else {
    DX <- (1/eta) * t(X)
    XDX <- X %*% DX
    M <- diag(N) + XDX/xi
    m <- solve(M, y)
    ymy <- t(y) %*% m
    if (s != 0) {
      new_M <- diag(N) + XDX/new_xi
      new_m <- solve(new_M, y)
      new_ymy <- t(y) %*% new_m
      cM <- diag(chol(M))^2
      new_cM <- diag(chol(new_M))^2
      curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w + ymy) -
        log((sqrt(xi) / psi) * (1 + xi * psi^2))
      new_ratio <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
        log((sqrt(new_xi) / psi) * (1 + new_xi * psi^2))
      acceptance_probability <- exp(new_ratio - curr_ratio + log(new_xi) -
                                      log(xi))
      u <- stats::runif(n = 1, min = 0, max = 1)
      if (u < acceptance_probability) {
        xi <- new_xi
        ymy <- new_ymy
        M <- new_M
      }
    }
  }
  # sigma update
  if(is.null(sigma2)){
    sigma2 <- 1
  } else{
    sigma2 <- 1/stats::rgamma(1, shape = (w + N)/2, rate = (w + ymy)/2)
  }
  # beta update
  diag_D <- 1/(eta * xi)
  u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
  f <- stats::rnorm(n = N, mean = 0, sd = 1)
  v <- X %*% u + f
  U <- diag_D * t(X)
  if (p < N) {
    yv <- (y)/sqrt(sigma2) - v
    Xyv <- t(X) %*% yv
    m <- solve(Q_star, Xyv)
    m_star <- yv - X %*% m
    new_beta <- sqrt(sigma2) * (u + U %*% m_star)
  } else {
    v_star <- solve(M, (y / sqrt(sigma2) - v))
    new_beta <- sqrt(sigma2) * (u + U %*% v_star)

  }
  eta <- Mhorseshoe:::rejection_sampler((new_beta^2)*xi/(2*sigma2), a, b)
  eta <- ifelse(eta <= 2.220446e-16, 2.220446e-16, eta)
  return(list(new_beta = new_beta, new_xi = xi, new_sigma2 = sigma2, new_eta = eta))
}
