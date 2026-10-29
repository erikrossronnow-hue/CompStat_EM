#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector e_step_t_cpp(NumericVector par, NumericVector y, double v) {
  int n = y.size();
  NumericVector w(n);
  double mu = par[0], s = v * par[1];
  for (int i = 0; i < n; ++i) {
    double r = y[i] - mu;
    w[i] = (v + 1) / (1 + r * r / s);
  }
  return w;
}
