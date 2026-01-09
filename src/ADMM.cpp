#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double soft_c(double a, double lambda)
{
  if (a > lambda)
  {
    return(a - lambda);
  }
  else if (a < -lambda)
  {
    return(a + lambda);
  }
  else
  {
    return(0);
  }
}

// [[Rcpp::export]]
Rcpp::List ADMM_ab_cpp(NumericVector C, NumericVector lambda, double mu=1, int maxit = 100, double eps = 1e-5) 
{
  int p = C.size();
  NumericVector b(p);
  NumericVector bo(p);
  double d = 0;
  double error = 1;
  int step = 0;
  double mu1 = 1/ (1+mu);
  
  while ( (error > eps) & (maxit > step)) 
  {
    for (int j = 0; j < p; j++) 
    {
      b[j] = soft_c(C[j] - mu*(sum(b)-b[j] + d), lambda[j])*mu1;
    }
    d += sum(b);
    error = mean(pow((b - bo), 2));
    step += 1;
    for (int i=0; i<p; i++) 
    {
      bo[i] = b[i]; 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("step") = step, Rcpp::Named("error") = error);
}

// [[Rcpp::export]]
Rcpp::List Soft_threshold_cpp(NumericVector C, NumericVector lambda) 
{
  int p = C.size();
  NumericVector b(p);
  for (int j = 0; j < p; j++) 
  {
    b[j] = soft_c(C[j], lambda[j]);
  }
  return Rcpp::List::create(Rcpp::Named("b") = b);
}

