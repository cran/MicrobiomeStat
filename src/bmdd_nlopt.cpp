// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#ifdef HAVE_NLOPT
#include <nlopt.hpp>

// Production-ready BMDD implementation using NLopt L-BFGS-B optimizer
// This provides the same convergence as nlminb with C++ speed

// Helper: compute alp0 and alp1 matrices from parameters (no covariates case)
inline void compute_alp_simple(const NumericVector& para, int m, int n,
                              NumericMatrix& alp0, NumericMatrix& alp1) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      alp0(i, j) = para[i];
      alp1(i, j) = para[m + i];
    }
  }
}

// Helper: compute pi matrix from parameters (no covariates case)
inline void compute_pi_simple(const NumericVector& para, int m, int n,
                             NumericMatrix& pi_mat) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      pi_mat(i, j) = para[i];
    }
  }
}

// E-step: Update gamma and beta
void e_step_count(NumericMatrix& gam, NumericMatrix& beta,
                 const NumericMatrix& alp0, const NumericMatrix& alp1,
                 const NumericMatrix& pi_mat, const NumericMatrix& W,
                 int m, int n, bool inner_loop, int inner_iterlim,
                 double inner_tol) {

  NumericMatrix alp_diff(m, n);
  NumericMatrix lg_alp_diff(m, n);
  NumericVector csum_alp_gam(n);

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      alp_diff(i, j) = alp0(i, j) - alp1(i, j);
      lg_alp_diff(i, j) = lgamma(alp0(i, j)) - lgamma(alp1(i, j));
    }
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double alp_gam_ij = alp0(i, j) - gam(i, j) * alp_diff(i, j);
      beta(i, j) = W(i, j) + alp_gam_ij;
    }
  }

  for (int j = 0; j < n; j++) {
    double sum_alp = 0;
    for (int i = 0; i < m; i++) {
      sum_alp += alp0(i, j) - gam(i, j) * alp_diff(i, j);
    }
    csum_alp_gam[j] = sum_alp;
  }

  NumericVector csum_beta(n);
  for (int j = 0; j < n; j++) {
    double sum = 0;
    for (int i = 0; i < m; i++) sum += beta(i, j);
    csum_beta[j] = sum;
  }

  if (inner_loop) {
    int iter = 1;

    while (iter <= inner_iterlim) {
      NumericMatrix gam_tmp = clone(gam);

      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          double tmp = csum_alp_gam[j] - alp0(i, j) + gam(i, j) * alp_diff(i, j);

          double h = lgamma(alp0(i, j) + tmp) - lgamma(alp1(i, j) + tmp) +
                     alp_diff(i, j) * (R::digamma(beta(i, j)) - R::digamma(csum_beta[j])) -
                     lg_alp_diff(i, j);

          gam(i, j) = 1.0 / (1.0 + (1.0 / pi_mat(i, j) - 1.0) * exp(h));

          double alp_gam_i = alp0(i, j) - gam(i, j) * alp_diff(i, j);
          csum_alp_gam[j] = tmp + alp_gam_i;

          csum_beta[j] -= beta(i, j);
          beta(i, j) = W(i, j) + alp_gam_i;
          csum_beta[j] += beta(i, j);
        }
      }

      double sum_sq_diff = 0.0, sum_sq_tmp = 0.0;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          double diff = gam(i, j) - gam_tmp(i, j);
          sum_sq_diff += diff * diff;
          sum_sq_tmp += gam_tmp(i, j) * gam_tmp(i, j);
        }
      }

      if (sum_sq_diff / sum_sq_tmp < inner_tol) {
        break;
      }

      iter++;
    }
  }
}

// Data structure for optimization
struct OptimData {
  NumericMatrix gam;
  NumericMatrix A;
  int m;
  int n;
};

// Objective function for NLopt with analytical gradients
double alpha_objective_nlopt(const std::vector<double>& x, std::vector<double>& grad, void* data) {
  OptimData* d = reinterpret_cast<OptimData*>(data);

  NumericVector para = wrap(x);
  NumericMatrix alp0(d->m, d->n), alp1(d->m, d->n);
  compute_alp_simple(para, d->m, d->n, alp0, alp1);

  double obj = 0.0;

  // Initialize gradient if needed
  if (!grad.empty()) {
    std::fill(grad.begin(), grad.end(), 0.0);
  }

  // Compute alp.gam = alp0 - gam * (alp0 - alp1) and column sums
  NumericMatrix alp_gam(d->m, d->n);
  NumericVector sum_alp_gam(d->n);

  for (int j = 0; j < d->n; j++) {
    double sum = 0.0;
    for (int i = 0; i < d->m; i++) {
      double alp_diff = alp0(i, j) - alp1(i, j);
      alp_gam(i, j) = alp0(i, j) - d->gam(i, j) * alp_diff;
      sum += alp_gam(i, j);
    }
    sum_alp_gam[j] = sum;
  }

  // Compute objective: sum(lgamma(colSums(alp.gam))) + sum(alp.gam * A - lg.alp.gam)
  for (int j = 0; j < d->n; j++) {
    obj += lgamma(sum_alp_gam[j]);

    for (int i = 0; i < d->m; i++) {
      double lg_alp_gam_ij = (1.0 - d->gam(i, j)) * lgamma(alp0(i, j)) +
                             d->gam(i, j) * lgamma(alp1(i, j));
      obj += alp_gam(i, j) * d->A(i, j) - lg_alp_gam_ij;
    }
  }

  // Compute gradient if needed
  if (!grad.empty()) {
    for (int j = 0; j < d->n; j++) {
      double digam_sum = R::digamma(sum_alp_gam[j]);

      for (int i = 0; i < d->m; i++) {
        // d/d(alp0[i]) = (1 - gam[i,j]) * (A[i,j] + digamma(sum) - digamma(alp0[i,j]))
        double grad_alp0 = (1.0 - d->gam(i, j)) * (d->A(i, j) + digam_sum - R::digamma(alp0(i, j)));

        // d/d(alp1[i]) = gam[i,j] * (A[i,j] + digamma(sum) - digamma(alp1[i,j]))
        double grad_alp1 = d->gam(i, j) * (d->A(i, j) + digam_sum - R::digamma(alp1(i, j)));

        grad[i] += -grad_alp0;  // Negative because NLopt minimizes
        grad[d->m + i] += -grad_alp1;
      }
    }
  }

  return -obj;  // Negative because NLopt minimizes
}

// [[Rcpp::export]]
List bmdd_nlopt(NumericMatrix W, String type,
               Nullable<NumericVector> para_alp_init = R_NilValue,
               Nullable<NumericVector> para_pi_init = R_NilValue,
               Nullable<NumericMatrix> gam_init = R_NilValue,
               int iterlim = 500, double tol = 1e-6, bool trace = false,
               bool inner_loop = true, int inner_iterlim = 20,
               double inner_tol = 1e-6,
               int alp_iterlim = 100, double alp_tol = 1e-6,
               double alp_min = 1e-3, double alp_max = 1e3) {

  int m = W.nrow();
  int n = W.ncol();

  // Check for NAs
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (NumericVector::is_na(W(i, j))) {
        stop("The OTU table contains NAs! Please remove!");
      }
    }
  }

  // Initialize parameters
  NumericMatrix gam(m, n);
  if (gam_init.isNotNull()) {
    gam = clone(as<NumericMatrix>(gam_init));
  } else {
    for (int i = 0; i < m * n; i++) {
      gam[i] = R::runif(0, 1);
    }
  }

  NumericVector para_alp(2 * m);
  if (para_alp_init.isNotNull()) {
    para_alp = clone(as<NumericVector>(para_alp_init));
  } else {
    for (int i = 0; i < m; i++) para_alp[i] = R::runif(0, 1);
    for (int i = m; i < 2 * m; i++) para_alp[i] = R::runif(1, 2);
  }

  // Ensure initial values satisfy bounds (NLopt requires this)
  for (int i = 0; i < 2 * m; i++) {
    if (para_alp[i] < alp_min) para_alp[i] = alp_min;
    if (para_alp[i] > alp_max) para_alp[i] = alp_max;
  }

  NumericVector para_pi(m);
  if (para_pi_init.isNotNull()) {
    para_pi = clone(as<NumericVector>(para_pi_init));
  } else {
    for (int i = 0; i < m; i++) {
      double sum = 0;
      for (int j = 0; j < n; j++) sum += gam(i, j);
      para_pi[i] = sum / n;
    }
  }

  NumericMatrix beta(m, n);
  NumericMatrix A(m, n);
  NumericMatrix alp0(m, n), alp1(m, n);
  NumericMatrix pi_mat(m, n);
  NumericMatrix W_work = clone(W);

  // Setup NLopt optimizer (will be reset for each M-step)
  // Using LBFGS without bounds first, will add bounds per iteration

  // Main EM loop
  int iter = 1;
  while (iter <= iterlim) {
    NumericVector para_alp_tmp = clone(para_alp);

    // Compute current alp0, alp1, pi
    compute_alp_simple(para_alp, m, n, alp0, alp1);
    compute_pi_simple(para_pi, m, n, pi_mat);

    // E-step
    e_step_count(gam, beta, alp0, alp1, pi_mat, W_work, m, n,
                 inner_loop, inner_iterlim, inner_tol);

    // Update A for M-step
    for (int j = 0; j < n; j++) {
      double csum = 0;
      for (int i = 0; i < m; i++) csum += beta(i, j);

      for (int i = 0; i < m; i++) {
        A(i, j) = R::digamma(beta(i, j)) - R::digamma(csum);
      }
    }

    // M-step: optimize alpha using NLopt
    OptimData optim_data;
    optim_data.gam = clone(gam);
    optim_data.A = clone(A);
    optim_data.m = m;
    optim_data.n = n;

    // Create new optimizer for this M-step
    nlopt::opt opt(nlopt::LD_LBFGS, 2 * m);

    std::vector<double> lower_bounds(2 * m, alp_min);
    std::vector<double> upper_bounds(2 * m, alp_max);
    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);

    opt.set_min_objective(alpha_objective_nlopt, &optim_data);
    // Match nlminb's default tolerances more closely
    opt.set_xtol_rel(alp_tol);     // Relative parameter tolerance
    opt.set_ftol_rel(alp_tol);     // Relative function tolerance
    opt.set_maxeval(alp_iterlim * 100);  // Allow many more evaluations

    std::vector<double> x(2 * m);
    for (int i = 0; i < 2 * m; i++) {
      // Clamp to bounds before optimization (safety check)
      x[i] = para_alp[i];
      if (x[i] < alp_min) x[i] = alp_min;
      if (x[i] > alp_max) x[i] = alp_max;
    }

    double minf;
    try {
      opt.optimize(x, minf);

      // Update para_alp
      for (int i = 0; i < 2 * m; i++) {
        para_alp[i] = x[i];
      }

    } catch(std::exception &e) {
      warning(std::string("NLopt optimization failed: ") + e.what());
    }

    // M-step: update pi
    for (int i = 0; i < m; i++) {
      double sum = 0;
      for (int j = 0; j < n; j++) sum += gam(i, j);
      para_pi[i] = sum / n;
    }

    if (trace) {
      Rcout << "Iteration " << iter << std::endl;
    }

    // Check convergence
    double sum_sq_diff = 0, sum_sq = 0;
    for (int i = 0; i < 2 * m; i++) {
      double diff = para_alp[i] - para_alp_tmp[i];
      sum_sq_diff += diff * diff;
      sum_sq += para_alp_tmp[i] * para_alp_tmp[i];
    }

    if (sum_sq_diff / sum_sq < tol) {
      break;
    }

    iter++;
  }

  // Final results
  compute_alp_simple(para_alp, m, n, alp0, alp1);
  compute_pi_simple(para_pi, m, n, pi_mat);

  // In no-covariate case, return alp0, alp1, and pi as vectors (matching R implementation)
  NumericVector alp0_vec(m);
  NumericVector alp1_vec(m);
  NumericVector pi_vec(m);

  for (int i = 0; i < m; i++) {
    alp0_vec[i] = para_alp[i];
    alp1_vec[i] = para_alp[m + i];
    pi_vec[i] = para_pi[i];
  }

  List alpha_list = List::create(Named("alp0") = alp0_vec, Named("alp1") = alp1_vec);

  return List::create(
    Named("gamma") = gam,
    Named("beta") = beta,
    Named("alpha") = alpha_list,
    Named("pi") = pi_vec,
    Named("para.alpha") = para_alp,
    Named("para.pi") = para_pi
  );
}

#else
// If NLopt is not available, provide a stub function that gives an informative error
// [[Rcpp::export]]
List bmdd_nlopt(NumericMatrix W, String type,
               Nullable<NumericVector> para_alp_init = R_NilValue,
               Nullable<NumericVector> para_pi_init = R_NilValue,
               Nullable<NumericMatrix> gam_init = R_NilValue,
               int iterlim = 500, double tol = 1e-6, bool trace = false,
               bool inner_loop = true, int inner_iterlim = 20,
               double inner_tol = 1e-6,
               int alp_iterlim = 100, double alp_tol = 1e-6,
               double alp_min = 1e-3, double alp_max = 1e3) {
  stop("NLopt library is not available. Please use bmdd() instead, or install NLopt and reinstall the package.");
  return List::create();  // Never reached
}
#endif
