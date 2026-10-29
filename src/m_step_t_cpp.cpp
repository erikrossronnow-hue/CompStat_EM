#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector m_step_t_cpp(NumericVector w, NumericVector y, double v) {
  int n = y.size();
  double sw = 0, sy = 0, s = 0;
  for (int i = 0; i < n; ++i) {
    sw += w[i];
    sy += w[i] * y[i];
  }
  double mu = sy / sw;
  for (int i = 0; i < n; ++i) {
    double r = y[i] - mu;
    s += w[i] * r * r;
  }
  return NumericVector::create(_["mu"] = mu, _["sigma2"] = s / (n * v));
}
