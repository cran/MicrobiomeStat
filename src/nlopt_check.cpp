#include <Rcpp.h>

// Simple function to check if NLopt compiled successfully
// [[Rcpp::export]]
bool test_nlopt_available() {
#ifdef HAVE_NLOPT
  return true;
#else
  return false;
#endif
}
