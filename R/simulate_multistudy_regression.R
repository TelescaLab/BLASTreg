#' Generate stacked multi-source regression data (target + sources)
#'
#' @description
#' Simulates a high-dimensional dataset for a **target** study followed by
#' \code{M} **source** studies, with a controllable subset of auxiliaries
#' designated as *informative*. The per-study coefficient vectors are built from
#' a shared sparse signal \eqn{\beta} and study-specific perturbations:
#' informative studies receive a shift of size \code{sig.delta1} on \eqn{h}
#' coordinates, while non-informative studies receive a shift of size
#' \code{sig.delta2} on \eqn{q} coordinates.
#'
#' Outcomes can be Gaussian (\eqn{y = X\beta + \varepsilon}) or logistic
#' (\eqn{y \sim \mathrm{Bernoulli}(\text{logit}^{-1}(X\beta))}).
#'
#' @details
#' - **Stacking/order:** Rows of \code{X} and entries of \code{y} are stacked as
#'   the target study **first**, then the \code{M} auxiliaries, in the order
#'   implied by \code{n.vec = c(n0, n1, ..., nM)}.
#' - **Base signal \eqn{\beta}:** first \code{s} entries set to \code{sig.beta},
#'   remaining \code{p - s} to 0.
#' - **Study perturbations:** Let \code{h = ceiling(p * contam_pct)}.
#'   Informative studies (first \code{size.A0}) shift \code{h} randomly chosen
#'   coordinates of \eqn{\beta} by \code{-sig.delta1}. Non-informative studies
#'   shift \code{q} coordinates by \code{-sig.delta2}, where
#'   \code{q = ceiling(ni_contam_factor * s)} if provided, otherwise \code{q = 2 * s}.
#' - **Design:** Each study draws rows independently from \eqn{N_p(0, I_p)}.
#' - **Test split:** If \code{test_set = TRUE}, a fraction \code{pct_test} of the
#'   **target** rows are held out and returned as \code{X.test}/\code{y.test};
#'   \code{n.vec[1]} is reduced accordingly.
#'
#' @param p Integer. Number of predictors.
#' @param s Integer. Number of non-zero signals in the base vector \eqn{\beta}.
#' @param M Integer. Number of source studies (total studies = \code{M + 1} including target).
#' @param size.A0 Integer in \code{[0, M]}. Number of *informative* source studies
#'   (by construction, these are the first \code{size.A0} auxiliaries).
#' @param n.vec Integer vector of length \code{M + 1} with per-study sample sizes
#'   \code{c(n0, n1, ..., nM)}; \code{n0} is the target sample size.
#' @param sig.beta Numeric. Magnitude of the non-zero entries in \eqn{\beta}.
#' @param sig.delta1 Numeric. Shift size applied to \code{h} coordinates for **informative** studies.
#' @param sig.delta2 Numeric. Shift size applied to \code{q} coordinates for **non-informative** studies.
#' @param contam_pct Numeric in \code{(0,1)}. Fraction of coordinates shifted in informative studies;
#'   used to set \code{h = ceiling(p * contam_pct)}. Default \code{0.015}.
#' @param ni_contam_factor Optional numeric. If supplied, sets
#'   \code{q = ceiling(ni_contam_factor * s)} for non-informative studies; otherwise \code{q = 2 * s}.
#' @param test_set Logical. If \code{TRUE}, hold out \code{pct_test} fraction of target rows as a test set.
#' @param pct_test Numeric in \code{(0,1)}. Fraction of target observations to hold out when \code{test_set = TRUE}.
#' @param type Character. Outcome model: \code{"gaussian"} (default) or \code{"logistic"}.
#'
#' @return
#' A list with components:
#'
#' - \code{X}: numeric matrix of size \code{sum(n.vec) × p} (or \code{(sum(n.vec) - n.test) × p} if a test split is made).
#' - \code{y}: numeric vector of length \code{sum(n.vec)} (or \code{sum(n.vec) - n.test} with test split).
#' - \code{coef.all}: list with
#'     - \code{beta0}: length-\code{p} base signal vector,
#'     - \code{W}: \code{p × M} matrix of source perturbations (one column per source study).
#' - \code{X.test}, \code{y.test}: held-out target design/response (or \code{NULL} if \code{test_set = FALSE}).
#' - \code{n.vec}: possibly updated per-study sizes (target reduced by \code{n.test} if a test split was made).
#'
#' @section Notes:
#' - Ensure \code{length(n.vec) == M + 1} and \code{size.A0 <= M}.
#' - This function draws from \code{mvtnorm::rmvnorm}; declare
#'   \verb{@importFrom mvtnorm rmvnorm} in your package if not already imported.
#' - Randomness is uncontrolled unless you set \code{set.seed()} beforehand.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' p <- 50; s <- 5; M <- 4
#' n.vec <- c(120, rep(100, M))
#'
#' sim <- simulate_multistudy_regression(
#'   p = p, s = s, M = M, size.A0 = 2, n.vec = n.vec,
#'   sig.beta = 0.7, sig.delta1 = 0.3, sig.delta2 = 1.0,
#'   contam_pct = 0.02, type = "gaussian"
#' )
#'
#' str(sim$X); length(sim$y)
#' str(sim$coef.all)
#'
#' # With a test split on the target study (~20% held out)
#' sim2 <- simulate_multistudy_regression(
#'   p = p, s = s, M = M, size.A0 = 2, n.vec = n.vec,
#'   sig.beta = 0.7, sig.delta1 = 0.3, sig.delta2 = 1.0,
#'   test_set = TRUE, pct_test = 0.2
#' )
#' c(train_target = sim2$n.vec[1], test_rows = nrow(sim2$X.test))
#' }
#'
#' @seealso [blast_select()], [blast_oracle()]
#' @importFrom mvtnorm rmvnorm
#' @export

simulate_multistudy_regression <- function(p, s, M, size.A0, n.vec, sig.beta, sig.delta1, sig.delta2,
                        contam_pct = 0.015, ni_contam_factor = NULL,  test_set = FALSE, pct_test = NULL,
                        type = "gaussian") {
  Sig.X <- diag(1, p)
  if(!is.null(ni_contam_factor)){
    q = ceiling(ni_contam_factor * s)
  } else{
    q = 2*s
  }
  h= ceiling(p * contam_pct)  #determines how many coefs we should bias for informative studies
  A0 = 1:size.A0  #this will make the first `size.A0` studies informative must be smaller than M
  beta0<-
    coef.all <-Coef.gen(s, h = h, q = q, size.A0 = size.A0,  M = M,   sig.beta = sig.beta,
                        sig.delta1 = sig.delta1, sig.delta2 = sig.delta2, p = p, exact=T)
  B <- cbind(coef.all$beta0, coef.all$W)
  beta0 <- coef.all$beta0
  #generate the data
  if(type == "gaussian"){
    X <- NULL
    y <- NULL
    for (k in 1:(M + 1)) {
      X <- rbind(X, mvtnorm::rmvnorm(n.vec[k], rep(0, p), Sig.X))
      ind.k <- ind.set(n.vec, k)
      y <- c(y, X[ind.k, ] %*% B[, k] + rnorm (n.vec[k], 0, 1))
    }

  } else if(type == "logistic"){
    X <- NULL
    y <- NULL
    for (k in 1:(M + 1)) {
      X <- rbind(X, mvtnorm::rmvnorm(n.vec[k], rep(0, p), Sig.X))
      ind.k <- ind.set(n.vec, k)
      y <- c(y, rbinom(n.vec[k], 1, 1/(1 + exp(-X[ind.k, ] %*% B[, k]))))
    }
  }

  X.test <- NULL
  y.test <- NULL
  n.test = NULL

  if(test_set){
    n0 <- n.vec[1]
    n.test = ceiling(n0 * pct_test)
    inds <- sample(1:n0, n.test)
    X.test <- X[inds, ]
    y.test <- y[inds]
    X = X[-inds, ]
    y = y[-inds]
    n.vec[1] = n.vec[1] - n.test
  }
  return(list(X = X, y = y, coef.all = coef.all, X.test = X.test, y.test = y.test, n.vec = n.vec))
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}


Coef.gen<- function(s, h,q=30, size.A0, M, sig.beta,sig.delta1, sig.delta2, p, exact=T){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- rep.col(beta0,  M)
  for(k in 1:M){
    if(k <= size.A0){
      if(exact){
        samp0<- sample(1:p, h, replace=F)
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1, h)
      }else{
        W[1:p,k] <-W[1:p,k] + rnorm(p, 0, h/p)
      }
    }else{
      if(exact){
        samp1 <- sample(1:p, q, replace = F)
        W[samp1,k] <- W[samp1,k] + rep(-sig.delta2,q)
      }else{
        W[1:p,k] <-W[1:p,k] + rnorm(p, 0, q/p)
      }
    }
  }
  return(list(W=W, beta0=beta0))
}
