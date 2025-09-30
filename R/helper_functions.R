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



