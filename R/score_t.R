score_t = function(par, y, nu) {
  r = y - par[1]
  s = par[2]
  d = nu * s + r^2
  cbind(mu = (nu + 1) * r / d, sigma2 = nu * (r^2 - s) / (2 * s * d))
}
