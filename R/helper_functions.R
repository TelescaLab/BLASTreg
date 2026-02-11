ind.set<- function(n.vec, k.vec){
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

set_data_inds <- function(gamma, n.vec){
  # Set up indices according to the gamma vector
  k.vec = which(gamma == 1) + 1
  ni.vec <- which(gamma == 0) + 1
  ind.kA <- ind.set(n.vec = n.vec, k.vec = k.vec)
  ni.vec <- which(gamma == 0) + 1
  ind.ni <- ind.set(n.vec = n.vec, k.vec = ni.vec)
  informative <- !is.null(ind.kA)
  noninformative <- !is.null(ind.ni)
  return(list(informative = informative, noninformative = noninformative,
              ind.kA = ind.kA, ind.ni = ind.ni))
}

normalize_log_probabilities <- function(log_p1, log_p0, scale = 1) {
  log_p1 <- scale * log_p1
  log_p0 <- scale * log_p0
  max_log_p <- max(log_p1, log_p0)
  exp_p1 <- exp(log_p1 - max_log_p)
  exp_p0 <- exp(log_p0 - max_log_p)
  norm_p1 <- exp_p1 / (exp_p1 + exp_p0)
  norm_p0 <- exp_p0 / (exp_p1 + exp_p0)
  return(c(norm_p1, norm_p0))
}

#' Precompute XtX and Xty for all source studies
precompute_XtXXty <- function(X, y, K, n.vec){
  list_XtX <- list()
  list_Xty <- list()
  for(k in 1:K){
    ind.k <- ind.set(n.vec = n.vec, k.vec = k + 1)
    X_k <- X[ind.k,]
    y_k <- y[ind.k]
    list_XtX[[k]] <- t(X_k) %*% X_k
    list_Xty[[k]] <- t(X_k) %*% y_k
  }
  return(list(list_XtX = list_XtX, list_Xty = list_Xty))
}

#' Construct candidate sets from sparsity index (Li et al., 2020)
construct_candidate_set <- function(X, y, n.vec) {
  M= length(n.vec)-1
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  ind.1<-ind.set(n.vec,1)
  for(k in 2: (M+1)){
    ind.k<-ind.set(n.vec,k)
    Xty.k <- t(X[ind.k,])%*%y[ind.k]/n.vec[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec[1]
    margin.T<-sort(abs(Xty.k),decreasing=T)
    Rhat[k] <-  sum(margin.T^2)
  }
  Tset <- which(Rhat[-1] <= quantile(Rhat[-1], probs = 0.25))
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use Rhat as the selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  return(Tset)
}


safe_log_det <- function(M, label = "") {
  M_sym <- 0.5 * (M + t(M))
  detM <- determinant(M_sym, logarithm = TRUE)
  # if (detM$sign <= 0 || !is.finite(detM$modulus)) {
  #   warning(sprintf("Matrix not positive definite in log_det(%s)", label))
  #   return(-1e10)
  # }
  if (!is.finite(detM$modulus)){
    dd <- -1e10
  }else{
    dd <- as.numeric(detM$modulus)
    }
  return(dd)
}

approx_marginal_likelihood <- function(
      y0, yA, yAb,
      X0, XA, XAb,
      d_A, d_Ab, d_delta, cut1=0.1)
{
  # constants and sample size
  log2pi <- log(2 * pi)
  n0     <- length(y0)
  nA     <- length(yA)
  nAb    <- length(yAb)

  # filtering effective zeros
  vp  <- (d_A + d_Ab)/2.0
  vp[sample(seq(1,ncol(X0)),3.0)] = 10.0
  #
  X0  <- X0[,vp > cut1]
  XA  <- XA[,vp > cut1]
  XAb <- XAb[,vp > cut1]
  #
  d_A <- d_A[vp > cut1]
  d_Ab <- d_Ab[vp > cut1]
  d_delta <-d_delta[vp > cut1]
  #
  p   <- ncol(X0)
  #
  # Loglik_Ab ------------------------------------------------
  M_Ab <- crossprod(XAb) + diag(1/d_Ab, p)
  M_Ab <- 0.5 * (M_Ab + t(M_Ab))
  Xty_Ab = t(XAb)%*%yAb
  tmp <- solve(M_Ab, Xty_Ab)
  quad_Ab <- drop(sum(yAb^2) - crossprod(Xty_Ab, tmp))
  #
  loglik_Ab <-
    lgamma(0.5 * (nAb + p - 1)) -
    0.5 * safe_log_det(M_Ab, "M_Ab") -
    0.5 * (nAb + p - 1) * log(0.5 * (quad_Ab + 1))
  # Loglik_A -------------------------------------------------
  Z0 = matrix(0, nrow=nA, ncol=p)
  Z  = rbind(cbind(X0, X0), cbind(XA, Z0))
  yy = c(y0, yA)
  M_Z <- crossprod(Z) + diag(1/c(d_A, d_delta), 2*p)
  M_Z <- 0.5 * (M_Z + t(M_Z))
  #
  Zty = t(Z)%*%yy
  #
  tmp <- solve(M_Z, Zty)
  quad_Z <- drop(sum(yy^2) - crossprod(Zty, tmp))
  #
  loglik_A <-
    lgamma(0.5 * (nA + n0 + 2*p - 1)) -
    0.5 * safe_log_det(M_Z, "M_Z") -
    0.5 * (nA + n0 + 2*p - 1) * log(0.5 * (quad_Z + 1))


  ## Marginal likelihood -------------------------------------
  approx_mlike = loglik_Ab + loglik_A

  return(approx_mlike)
}









