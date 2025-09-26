#' Bayesian Transfer Learning Shrinkage Regression with Source Selection
#'
#' @param X predictor matrix
#' @param y response vector
#' @param burn burn-in iterations
#' @param iter number of iterations
#' @param a parameter for rejection sampler
#' @param b parameter for rejection sampler
#' @param s tuning parameter for proposal distribution
#' @param tau global shrinkage parameter
#' @param sigma2 error variance
#' @param w inverse gamma prior shape parameter
#' @param alpha confidence level for credible intervals
#' @param n.vec vector of sample sizes for each study
#' @param iterEBstep number of iterations for empirical bayes step, if 0, no empirical bayes
#' @param sir whether to use SIR for empirical bayes
#'
#' @export
blast_select() <- function(X, y, n.vec, burn = 1000, iter = 3000, a = 1/5, b = 10,
                                  s = 0.8, tau = 1, sigma2 = 1, w = 1, alpha = 0.05,
                                  iterEBstep = 0, sir = TRUE, gamma_init = NULL, temp_scale = 1, adapt_temp_scale = FALSE){
  p <- ncol(X)
  N = sum(n.vec)
  K <- length(n.vec) - 1
  # initialize gamma using sparsity index
  if(is.null(gamma_init)){
    candidate_sets <- construct_candidate_set(X= X, y = y, n.vec = n.vec)
    gamma <- rep(0, K)
    gamma[candidate_sets[[1]]] <- 1
  } else{
    gamma <- gamma_init
  }

  # gamma <- gamma_init
  # Target Data setup
  n0 <- n.vec[1]
  N_k = N - n0
  X_0 <- X[1:n0, ]
  y_0 <- y[1:n0]
  eta_0 <- rep(1, p)
  xi_0 <- tau^(-2)
  Q0 <- t(X_0) %*% X_0
  sigma2_0 <- 1
  kappa_0 <- -0.5
  # Auxilliary Data setup
  etaA <- rep(1, p)
  xiA <- tau^(-2)
  eta_A_bar <- rep(1, p)
  xi_A_bar <- tau^(-2)
  wA_new <- rep(0, p)
  # gamma_prior_probs <- n.vec[2:(K+1)] / sum(n.vec[2:(K+1)])
  gamma_prior_probs <- rep(0.5, K)
  sigma2_A <- 1
  sigma2_A_bar <- 1
  kappa_A <- -0.5
  kappa_A_bar <- -0.5
  samp_size_A <- NULL
  samp_size_A_bar <- NULL
  # storage setup
  nmc <- burn + iter
  xi_0out <- rep(0, nmc)
  xiA_out <- rep(0, nmc)
  xiA_bar_out <- rep(0, nmc)
  sigma2_0out <- rep(0, nmc)
  sigma2_Aout <- rep(0, nmc)
  sigma2_A_bar_out <- rep(0, nmc)
  betaout <- matrix(0, nrow = nmc, ncol = p)
  wAout <- matrix(0, nrow = nmc, ncol = p)
  w0out <- matrix(0, nrow = nmc, ncol = p)
  deltaout <- matrix(0, nrow = nmc, ncol = p)
  gammaout <- matrix(0, nrow = nmc, ncol = K)
  etaout <- matrix(0, nrow = nmc, ncol = p)
  sigma2_0out <- rep(0, nmc)
  log_diff_values <- matrix(NA, nrow = nmc, ncol = K)
  eb_counter = 0
  # Alias for X0tX0
  X0_X0_t <- Q0
  X0_y0 <- t(X_0) %*% y_0
  # Precompute XtX Xty for aux studies
  precomp <- precompute_XtXXty(X, y, K, n.vec)
  EB <- iterEBstep != 0
  # Empirical bayes on the target data alone
  if(EB){
    EB_burn <- 1000
    EB_step_total <- iterEBstep + EB_burn
    # get horseshoe samples from target only
    eb_beta_samples <- matrix(0, nrow = EB_step_total + 1, ncol = p)
    eb_eta_samples <- matrix(0, nrow = EB_step_total + 1, ncol = p)
    eb_xi_samples <- rep(0, EB_step_total + 1)
    eb_sigma2_samples <- rep(0, EB_step_total + 1)
    # set initial values
    eb_beta_samples[1, ] <- rep(0, p)
    eb_eta_samples[1, ] <- rep(1, p)
    eb_xi_samples[1] <- xi_0
    eb_sigma2_samples[1] <- sigma2_0
    for(j in 2:(EB_step_total + 1)){
      targ_sample <- exact_horseshoe_gibbs_step(X_0, y_0, eta = eb_eta_samples[j-1, ], xi = eb_xi_samples[j-1],
                                                Q = Q0, w = w, s = s, a = a, b = b, p = p, sigma2 = eb_sigma2_samples[j-1],
                                                kappa = kappa_0, EB = FALSE)
      eb_beta_samples[j, ] <- targ_sample$new_beta
      eb_eta_samples[j, ] <- targ_sample$new_eta
      eb_xi_samples[j] <- targ_sample$new_xi
      eb_sigma2_samples[j] <- targ_sample$new_sigma2
    }
    tmp <- eb_em_max(kappa = kappa_0, N = n0, xi_out = eb_xi_samples[(EB_burn+1):(EB_step_total+1)],
                     sigma2_out = eb_sigma2_samples[(EB_burn+1):(EB_step_total+1)],
                     iterEBstep = iterEBstep, eb_counter = 0, sir = sir, i = iterEBstep)
    kappa_A <- tmp$kappa
    #kappa_A_bar <- tmp$kappa
    kappa_0 <- tmp$kappa
    if(tmp$kappa > 0){
      kappa_0 <- tmp$kappa / 10
      kappa_A_bar <- tmp$kappa *1.5
    } else{
      kappa_0 <- tmp$kappa * 1.2
      kappa_A_bar <- tmp$kappa / 2
    }

    print(kappa_A)
  }
  # Start MCMC
  for(i in 1:nmc){
    # Set up indices according to the gamma vector
    data_inds <- set_data_inds(gamma, n.vec)
    ind.kA <- data_inds$ind.kA
    ind.ni <- data_inds$ind.ni
    aux_num_ind.kA <- as.logical(gamma)
    aux_num_ind.ni <- !aux_num_ind.kA
    ### Individual updates
    # bias update
    delta_samp <- exact_horseshoe_gibbs_step(X_0, y_0 - X_0 %*% wA_new, eta = eta_0, xi = xi_0, Q = Q0,
                                             w = w, s = s, a = a, b = b, p = p, sigma2 = sigma2_0, kappa = kappa_0,
                                             EB = EB, kappa_0_trans = TRUE) #xiA_trunc = TRUE, xiA_trunc_val = xiA)
    #, trunc_val = log(xiA))
    # extract results for target
    delta_new <- delta_samp$new_beta
    xi_0 <- delta_samp$new_xi
    #xi_0 <- 50
    sigma2_0 <- delta_samp$new_sigma2
    # eta update target
    eta_0 <- delta_samp$new_eta
    ## sample w1, aux regression coefficient for informative samples
    if(data_inds$informative){
      X_A <- X[ind.kA, ]
      y_A <- y[ind.kA]
      newX <- rbind(X_0, X_A)
      newy <- c(y_0 - X_0 %*% delta_new, y_A)
      Q <- crossprod(newX)
      #Q_A <- add_matrices(precomp$list_XtX[aux_num_ind.kA])
    } else{
      # sample some random target data to fill in the gap
      pseudo_inds <- sample(1:n0, size = floor(0.05*n0))
      newX <- X[pseudo_inds, ]
      newy <- y[pseudo_inds]
      # newX <- rbind(X_0, X_A)
      # newy <- c(y_0 - X_0 %*% delta_new, y_A)
      # y_A = 0
      # X_A = t(rep(0, p))
      # newX <- t(rep(0, p))
      # newy <- 0
      # newX <- X_0
      # newy <- y_0 - X_0 %*% delta_new
      # Q <- crossprod(X_0)
      Q <- crossprod(newX)

    }
    # Auxiliary coefficient update for informative samples

    wA_samp <- exact_horseshoe_gibbs_step(X = newX, y = newy, eta = etaA, xi = xiA, Q = Q,
                                          w = w, s = s, a = a, b = b, p = p, sigma2 = sigma2_A, kappa = kappa_A,
                                          EB = EB)
    # extract results for informative
    wA_new <- wA_samp$new_beta
    xiA <- wA_samp$new_xi
    #xiA <- 100
    # xiA <-sum(n.vec[2:(K+1)]) / 4 # fix at "truth" -- tau = no. of signals / sample size
    sigma2_A <- wA_samp$new_sigma2
    # eta update informative
    etaA <- wA_samp$new_eta
    #eta_A <- c(rep(0.1, 16), rep(5, p-16))

    wAout[i, ] <- wA_new

    ## sample w0, aux regression coefficient for non-informative samples
    if(data_inds$noninformative){
      # if non-informative set is not empty, then access the data at the indices
      X_A_bar <- X[ind.ni, ]
      y_A_bar <- y[ind.ni]
      Q_A_bar <- t(X_A_bar) %*% X_A_bar
      #Q_A_bar <- add_matrices(precomp$list_XtX[aux_num_ind.ni])
    } else{
      # if non-informative set is empty, set y_A_bar and X_A_bar to 0 vectors
      # randomly sample some auxiliary data to fill in the gap
      pseudo_inds <- sample((n0 + 1):N, size = floor(0.03*N_k))
      y_A_bar <- y[sample(pseudo_inds)] # shuffled y's to make it more realistic
      X_A_bar <- X[pseudo_inds, ]
      # y_A_bar = 0
      # X_A_bar = t(rep(0, p))
      Q_A_bar <- t(X_A_bar) %*% X_A_bar
    }

    w0_samp <- exact_horseshoe_gibbs_step(X_A_bar, y_A_bar, eta = eta_A_bar, xi = xi_A_bar, Q = Q_A_bar,
                                          w = w, s = s, a = a, b = b, p = p, sigma2 = sigma2_A_bar, kappa = kappa_A_bar,
                                          EB = EB, kappa_A_bar_trans = TRUE)
    # extract results for non-informative
    w0_new <- w0_samp$new_beta
    xi_A_bar <- w0_samp$new_xi
    sigma2_A_bar <- w0_samp$new_sigma2
    # eta update non-informative
    eta_A_bar <- w0_samp$new_eta
    w0out[i, ] <- w0_new

    deltaout[i, ] <- delta_new
    # Sample beta
    betaout[i, ] <- wA_new + delta_new
    # Prepare prior covariance matrices for posterior gammma probability calculation
    # Sigma_A <- sigma2_A * (1 / xiA) * diag(1 / etaA)
    # Sigma_delta <- sigma2_0 * (1 / xi_0) * diag(1 / eta_0)
    # Sigma_A_bar <- sigma2_A_bar * (1 / xi_A_bar) * diag(1 / eta_A_bar)
    Sigma_A <-  (1 / xiA) * diag(1 / etaA)
    Sigma_delta <- (1 / xi_0) * diag(1 / eta_0)
    Sigma_A_bar <- (1 / xi_A_bar) * diag(1 / eta_A_bar)
    # Precompute XtX and Xty's
    # Compute the posterior probability of gamma for k = 1, ..., K using marginal likelihood method
    # update order for gamma
    #update_order <- sample(1:K)
    #update_order <- 1:K
    for(k in 1:K){
      gamma_tmp <- gamma
      # Set appropriate indices where gamma_k = 1
      gamma_tmp[k] <- 1
      inds.tmp <- set_data_inds(gamma_tmp, n.vec)
      # Check whether the non-informative set is empty
      # If it is, set y_A_bar and X_A_bar to 0
      # Note that informative set is not empty in this case
      if(!inds.tmp$noninformative){
        y_A_bar = 0
        X_A_bar = t(rep(0, p))
        # pseudo_inds <- sample((n0 + 1):N, size = floor(0.02*N_k))
        # y_A_bar <- y[sample(pseudo_inds)] # shuffled y's
        # X_A_bar <- X[pseudo_inds, ]
      } else{
        y_A_bar = y[inds.tmp$ind.ni]
        X_A_bar = X[inds.tmp$ind.ni, ]
      }
      aux_num_ind.kA <- as.logical(gamma_tmp)
      aux_num_ind.ni <- !aux_num_ind.kA
      precomp_XtX_A <- precomp$list_XtX[aux_num_ind.kA]
      precomp_Xty_A <- precomp$list_Xty[aux_num_ind.kA]
      precomp_XtX_A_bar <- precomp$list_XtX[aux_num_ind.ni]
      precomp_Xty_A_bar <- precomp$list_Xty[aux_num_ind.ni]
      gamma_prior_prob = gamma_prior_probs[k]

      log_p1 <- log_marginal_likelihood_marg_sig(
        y0 = y_0,
        yA = y[inds.tmp$ind.kA],
        yAb = y_A_bar,
        X0 = X_0,
        XA = X[inds.tmp$ind.kA, ],
        XAb = X_A_bar,
        #delta = delta_new,
        d_A = diag(Sigma_A),
        d_Ab = diag(Sigma_A_bar),
        d_delta = diag(Sigma_delta),
        X0_X0_t = X0_X0_t,
        X0_y0 = X0_y0,
        precomp_XtX_A = precomp_XtX_A,
        precomp_Xty_A = precomp_Xty_A,
        precomp_XtX_Ab = precomp_XtX_A_bar,
        precomp_Xty_Ab = precomp_Xty_A_bar,
        informative = inds.tmp$informative,
        noninformative = inds.tmp$noninformative
      )
      # Set appropriate indices where gamma_k = 0
      gamma_tmp[k] <- 0
      aux_num_ind.kA <- as.logical(gamma_tmp)
      aux_num_ind.ni <- !aux_num_ind.kA
      precomp_XtX_A <- precomp$list_XtX[aux_num_ind.kA]
      precomp_Xty_A <- precomp$list_Xty[aux_num_ind.kA]
      precomp_XtX_A_bar <- precomp$list_XtX[aux_num_ind.ni]
      precomp_Xty_A_bar <- precomp$list_Xty[aux_num_ind.ni]
      inds.tmp <- set_data_inds(gamma_tmp, n.vec)
      # Follow the same procedure as above
      if(!inds.tmp$informative){
        y_A = 0
        X_A = t(rep(0, p))
      } else{
        y_A = y[inds.tmp$ind.kA]
        X_A = X[inds.tmp$ind.kA, ]
      }
      log_p0 <- log_marginal_likelihood_marg_sig(
        y0 = y_0,
        yA = y_A,
        yAb = y[inds.tmp$ind.ni],
        X0 = X_0,
        XA = X_A,
        XAb = X[inds.tmp$ind.ni, ],
        #delta = delta_new,
        d_A = diag(Sigma_A),
        d_Ab = diag(Sigma_A_bar),
        d_delta = diag(Sigma_delta),
        X0_X0_t = X0_X0_t,
        X0_y0 = X0_y0,
        precomp_XtX_A = precomp_XtX_A,
        precomp_Xty_A = precomp_Xty_A,
        precomp_XtX_Ab = precomp_XtX_A_bar,
        precomp_Xty_Ab = precomp_Xty_A_bar,
        informative = inds.tmp$informative,
        noninformative = inds.tmp$noninformative
      )
      # Calculate the acceptance probability
      # Normalize the log probabilities for numerical stability
      probs <- normalize_log_probabilities(log_p1, log_p0, scale = temp_scale)
      p1 <- probs[1]
      p0 <- probs[2]
      gamma_old <- gamma[k]
      gamma[k] <- rbinom(1, 1, p1)
      # if(all(gamma == 0)){
      #   gamma[k] <- 1
      # }
      # if(all(gamma == 1)){
      #   gamma[k] <- 0
      # }
      if(i %% 100 == 0){
        cat("Gamma: ", k, "\n")
        cat("Log P1: ", log_p1, "\n")
        cat("Log P0: ", log_p0, "\n")
        cat("Diff: ", log_p1 - log_p0, "\n")
        cat("P1: ", p1, "\n")
      }
    }
    # save results
    gammaout[i, ] <- gamma
    etaout[i, ] <- eta_0
    xi_0out[i] <- xi_0
    xiA_out[i] <- xiA
    xiA_bar_out[i] <- xi_A_bar
    sigma2_0out[i] <- sigma2_0
    sigma2_Aout[i] <- sigma2_A
    sigma2_A_bar_out[i] <- sigma2_A_bar


    # print every 100 iterations
    if(i %% 100 == 0){
      cat("Iteration: ", i, "\n")
    }
  }
  betaout <- betaout[(burn+1):nrow(betaout), ]
  gammaout <- gammaout[(burn+1):nmc, ]
  lambdaout <- 1/sqrt(etaout[(burn+1):nmc, ])
  tauout <- 1/sqrt(xi_0out[(burn+1):nmc])
  sigma2_0out <- sigma2_0out[(burn+1):nmc]
  betahat <- apply(betaout, 2, mean)
  betahat_median <- apply(betaout, 2, median)
  lambdahat <- apply(lambdaout, 2, mean)
  tauhat <- mean(tauout)
  sigma2hat <- mean(sigma2_0out)
  gammasums <- colSums(gammaout)
  leftci <- apply(betaout, 2, stats::quantile, probs = alpha/2)
  rightci <- apply(betaout, 2, stats::quantile, probs = 1-alpha/2)
  result <- list(BetaHat = betahat, BetaHatMedian = betahat_median, LeftCI = leftci, RightCI = rightci,
                 Sigma2Hat = sigma2hat, TauHat = tauhat, LambdaHat = lambdahat,
                 BetaSamples = betaout, LambdaSamples = lambdaout,
                 TauSamples = tauout, Sigma2Samples = sigma2_0out,
                 GammaSamples = gammaout, W0Samples = w0out, WA_Samples = wAout,
                 GammaSums = gammasums)
  return(result)

}
