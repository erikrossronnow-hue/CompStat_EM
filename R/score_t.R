score_t = function(par, y, v) {
  r = y - par[1]
  s = par[2]
  d = v * s + r^2
  cbind(mu = (v + 1) * r / d, sigma2 = v * (r^2 - s) / (2 * s * d))
}
