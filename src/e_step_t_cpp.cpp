#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector e_step_t_cpp(NumericVector par, NumericVector y, double nu) {
  int n = y.size();
  NumericVector w(n);
  double mu = par[0], s = nu * par[1];
  for (int i = 0; i < n; ++i) {
    double r = y[i] - mu;
    w[i] = (nu + 1) / (1 + r * r / s);
  }
  return w;
}
