#' Generates data for high-dimensional linear regression
#'
#' @param p number of predictors
#' @param s number of signals
#' @param M number of studies
#' @param size.A0 number of informative studies
#' @param n.vec vector of sample sizes for each study
#' @param sig.beta signal size
#' @param beta_bias signal bias
#' @param contam_pct contamination percentage
#' @param ni_contam_factor non-informative contamination factor
#' @param test_set whether to generate a test set
#' @param pct_test percentage of data to use for test set
#' @param type type of data to generate
#'
#' @export
simulate_multistudy_regression <- function(p, s, M, size.A0, n.vec, sig.beta, sig.delta1, sig.delta2, beta_bias,
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
  #W[1,]<-W[1,]-2*sig.beta
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
