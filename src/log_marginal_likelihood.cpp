#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Robust log(det) computation with symmetry enforcement
double safe_log_det2(const arma::mat& M, const std::string& label) {
  arma::mat M_sym = 0.5 * (M + M.t());
  double val, sign;
  arma::log_det(val, sign, M_sym);
  // if (sign <= 0 || !std::isfinite(val)) {
  //   Rcpp::Rcout << "Warning: Matrix not positive definite in log_det(" << label << ")\n";
  //   return -1e10;
  // }
  return val;
}

// Helper to sum list of matrices
arma::mat add_precomp_matrices2(const List& mats) {
  arma::mat sum = as<arma::mat>(mats[0]);
  for (int i = 1; i < mats.size(); ++i) {
    sum += as<arma::mat>(mats[i]);
  }
  return sum;
}

arma::vec add_vectors2(const List& vecs) {
  arma::vec sum = as<arma::vec>(vecs[0]);
  for (int i = 1; i < vecs.size(); ++i) {
    sum += as<arma::vec>(vecs[i]);
  }
  return sum;
}

// [[Rcpp::export]]
double log_marginal_likelihood_marg_sig(
    const arma::vec& y0,
    const arma::vec& yA,
    const arma::vec& yAb,
    const arma::mat& X0,
    const arma::mat& XA,
    const arma::mat& XAb,
    const arma::vec& d_A,
    const arma::vec& d_Ab,
    const arma::vec& d_delta,
    const arma::mat& X0_X0_t,
    const arma::vec& X0_y0,
    const List& precomp_XtX_A,
    const List& precomp_Xty_A,
    const List& precomp_XtX_Ab,
    const List& precomp_Xty_Ab,
    bool informative,
    bool noninformative
) {
  const double log2pi = std::log(2.0 * M_PI);
  int n0 = y0.n_elem, nA = yA.n_elem, nAb = yAb.n_elem;
  int p = X0.n_cols;
  const double jitter = 1e-8;

  arma::mat XtX_A, XtX_Ab;
  arma::vec Xty_A, Xty_Ab;

  if (!informative) {
    XtX_A = XA.t() * XA;
    Xty_A = XA.t() * yA;
  } else {
    XtX_A = add_precomp_matrices2(precomp_XtX_A);
    Xty_A = add_vectors2(precomp_Xty_A);
  }

  if (!noninformative) {
    XtX_Ab = XAb.t() * XAb;
    Xty_Ab = XAb.t() * yAb;
  } else {
    XtX_Ab = add_precomp_matrices2(precomp_XtX_Ab);
    Xty_Ab = add_vectors2(precomp_Xty_Ab);
  }

  arma::vec inv_dA = 1.0 / (d_A + jitter);
  arma::vec inv_dAb = 1.0 / (d_Ab + jitter);
  arma::vec inv_dDelta = 1.0 / (d_delta + jitter);



  arma::mat M_Ab = XtX_Ab + diagmat(inv_dAb);
  M_Ab = 0.5 * (M_Ab + M_Ab.t());
  double quad_Ab = dot(yAb, yAb) - as_scalar(Xty_Ab.t() * solve(M_Ab, Xty_Ab));
  double loglik_Ab =
    -0.5 * nAb * log2pi
    //-0.5 * nAb * std::log(sigma2_Ab)
    //-0.5 * p * std::log(sigma2_Ab)
    + lgamma(0.5 * (nAb + p - 1))
    -0.5 * sum(log(d_Ab + jitter))
    -0.5 * safe_log_det2(M_Ab, "M_Ab")
    -0.5 * (nAb + p -1) * log(0.5 * (quad_Ab + 1));


    arma::mat M_A = X0_X0_t + XtX_A + diagmat(inv_dA);
    M_A = 0.5 * (M_A + M_A.t());
    arma::mat M_A_X0_X0_t = solve(M_A, X0_X0_t);
    arma::mat M_delta = X0_X0_t +  diagmat(inv_dDelta) - X0_X0_t * M_A_X0_X0_t;
    M_delta = 0.5 * (M_delta + M_delta.t());
    arma::mat M_delta_inv = inv(M_delta);



    double term1 =
      -0.5 * (nA) * log2pi
      //-0.5 * dot(yA, yA) / sigma2_A
      //-0.5 * p * std::log(sigma2_A)
      //-0.5 * p * std::log(sigma2_0)
      -0.5 * sum(log(d_A + jitter))
      -0.5 * sum(log(d_delta + jitter))
      -0.5 * safe_log_det2(M_A, "M_A")
      -0.5 * safe_log_det2(M_delta, "M_delta");

      arma::vec M_A_X0_y0 = solve(M_A, X0_y0);
      arma::vec M_A_XA_yA = solve(M_A, Xty_A);
      double quad_A = dot(yA, yA) - as_scalar(Xty_A.t() * M_A_XA_yA)
        + dot(y0, y0) - as_scalar(X0_y0.t() * M_A_X0_y0) - 2 * as_scalar(Xty_A.t() * M_A_X0_y0);
      arma::vec b0 = X0_y0 - X0_X0_t * M_A_XA_yA - X0_X0_t * M_A_X0_y0;
      double quad_A_bMb = quad_A - as_scalar(b0.t() * M_delta_inv * b0);
      double loglik_A = lgamma(0.5 * (n0 + nA + 4*p -1)) - 0.5 * (n0 + nA + 4*p -1) * log(0.5 * (quad_A_bMb + 1));


      // double term2 = 0.5 / (sigma2_A * sigma2_A) * as_scalar(Xty_A.t() * solve(M_A, Xty_A));
      // double term3 = -0.5 * dot(y0, y0) + 0.5 / (sigma2_0 * sigma2_0) * as_scalar(X0_y0.t() * solve(M_A, X0_y0));
      // double term4 = 1 / (sigma2_A * sigma2_0) * as_scalar(Xty_A.t() * solve(M_A, X0_y0));





      //arma::mat temp = solve(M_A, diff);
      //arma::mat inner = XtX_0_scaled * solve(M_delta, XtX_0_scaled * temp);
      // arma::rowvec temp = Xty_sum_scaled.t() * M_A_X0_X0_t;
      // double quad_term = as_scalar(temp * M_delta_inv * temp.t());
      // double term5 = 0.5 * as_scalar(Xty_0_scaled.t() * M_delta_inv * Xty_0_scaled) - 0.5 * quad_term;
      // double logdet_Mdelta = safe_log_det1(M_delta, "M_delta");

      return loglik_Ab + term1 + loglik_A;
}
